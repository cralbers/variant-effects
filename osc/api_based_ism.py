import os
import csv
import time
import threading
import multiprocessing
import queue
import requests
import numpy as np
from pathlib import Path
from dotenv import load_dotenv

load_dotenv(dotenv_path=Path(__file__).parent / ".env")

GENE_STRAND = "+"
GENE_NAME   = "LMNA"
CHROM       = "chr1"
ISM_START = 156_082_572    
ISM_END   = 156_082_575
CHUNK_BP        = 500   # bp per score_ism_variants() call
SEQ_LENGTH      = 1_048_576 #model input/context window
CALL_TIMEOUT_S  = 1000
MAX_WORKERS     = 16    # parallelism on AG's end. internal score_ism_variants concurrency per call
TISSUE_ONTOLOGY = "UBERON:0000948"  # heart
SCORES_CSV = str(Path(__file__).parent / f"bag3_ism_{ISM_START}_{ISM_END}_{SEQ_LENGTH}.csv")
# Dynamically picks up every non-empty AG_API_KEY_<N> env var present, in
# numeric order -- add or remove keys in .env and the worker count follows,
API_KEYS = [v for _, v in sorted(
    ((k, os.environ[k].strip()) for k in os.environ
     if k.startswith("AG_API_KEY_") and k.removeprefix("AG_API_KEY_").isdigit()),
    key=lambda kv: int(kv[0].removeprefix("AG_API_KEY_")),
) if v]

#Note: To add a scorer/modality, _build_scorers() and extract_score() will need updating

# RECOMMENDED columns are pulled straight from RECOMMENDED_VARIANT_SCORERS;
# everything else is COMPLEMENT, even if it happens to match the recommended one.
def _build_scorers():
    from alphagenome.models import dna_client, variant_scorers

    return [
        ("SPLICE_JUNCTIONS, RECOMMENDED", "SPLICE_JUNCTIONS",
         variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"]),
        ("SPLICE_SITES, RECOMMENDED", "SPLICE_SITES",
         variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"]),
        ("SPLICE_SITES, COMPLEMENT", "SPLICE_SITES",
         variant_scorers.CenterMaskScorer(
             requested_output=dna_client.OutputType.SPLICE_SITES, width=501,
             aggregation_type=variant_scorers.AggregationType.DIFF_SUM)),
        ("SPLICE_SITE_USAGE, RECOMMENDED", "SPLICE_SITE_USAGE",
         variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"]),
        ("SPLICE_SITE_USAGE, COMPLEMENT", "SPLICE_SITE_USAGE",
         variant_scorers.CenterMaskScorer(
             requested_output=dna_client.OutputType.SPLICE_SITE_USAGE, width=501,
             aggregation_type=variant_scorers.AggregationType.DIFF_LOG2_SUM)),
    ]

def fetch_reference_sequence(chrom_bare: str, start_1based: int, end_1based: int) -> str:
    url = (f"https://rest.ensembl.org/sequence/region/human/"
           f"{chrom_bare}:{start_1based}-{end_1based}")
    r = requests.get(url, headers={"Content-Type": "application/json"}, timeout=30)
    r.raise_for_status()
    return r.json()["seq"].upper()

def build_all_variant_rows(ref_seq: str, start_1based: int, score_cols: list) -> list:
    blank_scores = {c: "" for c in score_cols}
    return [
        {"position": start_1based + i, "ref": ref, "alt": alt, **blank_scores}
        for i, ref in enumerate(ref_seq)
        for alt in ([b for b in "ACGT" if b != ref] if ref in "ACGT" else list("ACGT"))
    ]

# extract_score's shape guards double as emptiness guards: an empty adata.X
# yields an all-False col_mask/row_mask below, which the .any() checks catch.
def extract_score(adata, output_type_name: str) -> float:
    X, var = adata.X, adata.var
    if output_type_name == "SPLICE_SITES":
        col_mask = ((var["name"] == "donor") | (var["name"] == "acceptor")) & \
                   (var["strand"] == GENE_STRAND)
    elif output_type_name == "SPLICE_JUNCTIONS":
        col_mask = (var["ontology_curie"] == TISSUE_ONTOLOGY) & \
                   (var["Assay title"] == "total RNA-seq")
    else:  # SPLICE_SITE_USAGE
        col_mask = (var["ontology_curie"] == TISSUE_ONTOLOGY) & \
                   (var["Assay title"] == "total RNA-seq") & \
                   (var["strand"] == GENE_STRAND)
    col_mask = col_mask.values
    if not col_mask.any():
        return float("nan")

    obs = adata.obs
    row_mask = (obs["gene_name"] == GENE_NAME).values if "gene_name" in obs.columns \
        else np.ones(X.shape[0], dtype=bool)
    if not row_mask.any():
        return float("nan")

    sub = np.asarray(X)[np.ix_(row_mask, col_mask)].flatten().astype(float)
    finite = sub[np.isfinite(sub)]
    if finite.size == 0:
        return float("nan")
    return float(finite[np.argmax(np.abs(finite))])

def run_one_piece(api_key, chunk_start, chunk_end, out_queue):
    from alphagenome.data import genome
    from alphagenome.models import dna_client

    scorers = _build_scorers()
    model = dna_client.create(api_key)
    ism_interval = genome.Interval(CHROM, chunk_start, chunk_end)
    seq_interval = ism_interval.resize(SEQ_LENGTH)

    t0 = time.time()
    try:
        result = model.score_ism_variants(
            interval=seq_interval,
            ism_interval=ism_interval,
            variant_scorers=[s for _, _, s in scorers],
            organism=dna_client.Organism.HOMO_SAPIENS,
            progress_bar=False,
            max_workers=MAX_WORKERS,
        )
        partial = {}
        for variant_scores in result:
            v = variant_scores[0].uns["variant"]
            key = (v.position, v.reference_bases, v.alternate_bases)
            row = {}
            for (col_name, ot_name, _scorer), adata in zip(scorers, variant_scores):
                score = extract_score(adata, ot_name)
                row[col_name] = "NaN" if np.isnan(score) else score
            partial[key] = row
        out_queue.put({"status": "ok", "elapsed": time.time() - t0, "partial": partial})
    except Exception:
        out_queue.put({"status": "error", "elapsed": time.time() - t0, "partial": {}})

def key_dispatcher(worker_id, api_key, my_chunks, rows_by_position,
                    lock, write_csv_fn, counters):
    ctx = multiprocessing.get_context("spawn")

    for chunk_start, chunk_end in my_chunks:
        out_queue = ctx.Queue()
        p = ctx.Process(target=run_one_piece,
                         args=(api_key, chunk_start, chunk_end, out_queue),
                         daemon=True)
        p.start()
        try:
            result = out_queue.get(timeout=CALL_TIMEOUT_S)
        except queue.Empty:
            result = {"status": "timeout", "elapsed": CALL_TIMEOUT_S, "partial": {}}
        p.join(timeout=5)
        if p.is_alive():
            p.terminate()

        ok = result["status"] == "ok"
        with lock:
            if ok:
                for (pos, ref, alt), row in result["partial"].items():
                    for r in rows_by_position.get(pos, []):
                        if r["ref"] == ref and r["alt"] == alt:
                            r.update(row)
                write_csv_fn()
            counters["ok" if ok else "failed"] += 1
            done = counters["ok"] + counters["failed"]

        status_msg = f"ok  {result['elapsed']:.1f}s" if ok else \
            "FAILED -- a rerun will likely fix it"
        print(f"  [key {worker_id}] {CHROM}:{chunk_start}-{chunk_end}  "
              f"{status_msg}  ({done}/{counters['total']})", flush=True)

def _write_csv(rows_list, score_cols):
    tmp = SCORES_CSV + ".tmp"
    with open(tmp, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["position", "ref", "alt"] + score_cols)
        w.writeheader()
        w.writerows(rows_list)
    os.replace(tmp, SCORES_CSV)

def main() -> None:
    scorers = _build_scorers()
    score_cols = [c for c, _, _ in scorers]

    if not Path(SCORES_CSV).exists():
        print(f"Fetching reference sequence ({CHROM}:{ISM_START+1}-{ISM_END})...")
        ref_seq = fetch_reference_sequence(CHROM.replace("chr", ""), ISM_START + 1, ISM_END)
        rows_list = build_all_variant_rows(ref_seq, ISM_START + 1, score_cols)
        _write_csv(rows_list, score_cols)
        print(f"Pre-initialised {SCORES_CSV} with {len(rows_list):,} rows "
              f"({len(ref_seq):,} positions x 3 alts).")
    else:
        with open(SCORES_CSV, newline="") as f:
            reader = csv.DictReader(f)
            missing_cols = [c for c in score_cols if c not in reader.fieldnames]
            if missing_cols:
                raise SystemExit(
                    f"{SCORES_CSV} doesn't have column(s) {missing_cols} -- its "
                    f"scorer config doesn't match the current _build_scorers(). "
                    f"Delete the file or change the region to start fresh."
                )
            rows_list = [{**r, "position": int(r["position"])} for r in reader]
        print(f"Loaded existing {SCORES_CSV} ({len(rows_list):,} rows).")

    rows_by_position: dict = {}
    for r in rows_list:
        rows_by_position.setdefault(r["position"], []).append(r)

    def chunk_is_done(cs, ce):
        return not any(r[c] == ""
                        for pos in range(cs + 1, ce + 1)
                        for r in rows_by_position.get(pos, [])
                        for c in score_cols)

    chunks = [(p, min(p + CHUNK_BP, ISM_END)) for p in range(ISM_START, ISM_END, CHUNK_BP)]
    todo = [c for c in chunks if not chunk_is_done(*c)]
    print(f"{len(chunks)} chunks total | {len(chunks) - len(todo)} already done | "
          f"{len(todo)} to run\n")
    if not todo:
        print("All chunks complete.")
        return

    lock = threading.Lock()
    write_csv_fn = lambda: _write_csv(rows_list, score_cols)
    key_chunks = {i: todo[i::len(API_KEYS)] for i in range(len(API_KEYS))}
    counters = {"ok": 0, "failed": 0, "total": len(todo)}

    print(f"Launching {len(API_KEYS)} key-dispatchers over {len(todo)} chunks "
          f"({CHUNK_BP} bp each)...\n")
    run_start = time.time()

    threads = [
        threading.Thread(
            target=key_dispatcher,
            args=(i + 1, key, key_chunks[i], rows_by_position, lock, write_csv_fn, counters),
            daemon=True,
        )
        for i, key in enumerate(API_KEYS)
    ]
    for t in threads:
        t.start()
    for t in threads:
        t.join()

    elapsed_min = (time.time() - run_start) / 60
    print()
    print("=" * 60)
    print(f"Finished in {elapsed_min:.1f} min | {counters['ok']} chunks ok | "
          f"{counters['failed']} chunks left blank")
    print("Re-run the script to retry the blank cells." if counters["failed"]
          else "All chunks complete.")


if __name__ == "__main__":
    main()
