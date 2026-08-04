from alphagenome_research.model import dna_model
from alphagenome import colab_utils
from alphagenome.data import gene_annotation
from alphagenome.data import genome
from alphagenome.data import transcript as transcript_utils
from alphagenome.interpretation import ism
from alphagenome.models import dna_client
from alphagenome.models import variant_scorers
from alphagenome.visualization import plot_components
import matplotlib.pyplot as plt
import pandas as pd
import sys
import jax
from numba import cuda
import pickle
import matplotlib.patches as mpatches
import numpy as np
from pysam import VariantFile
from tqdm import tqdm
import os
import gc
import math

if len(sys.argv) != 3:
    print("Usage: python run_ism_regulatory_regions.py <bed_file_path> <output_folder>")
    print("")
    sys.exit(1)

LMNA_START = 156_082_572
LMNA_END = 156_140_081
gene_symbol = "LMNA"
LMNA_INTERVAL = genome.Interval('chr1', 156_082_572, 156_140_081)


BASE_PATH = '/users/PAS2905/coraalbers/'
AG_DATA_PATH = '/users/PAS2905/coraalbers/ag/ag_data/'

HG38_FASTA_PATH = '/users/PAS2905/coraalbers/ag/hg38.fa'
HG38_GTF_PATH = '/users/PAS2905/coraalbers/ag/ag_data/gencode.v46.annotation.gtf.gz.feather'
HG38_SPLICE_START_PATH = '/users/PAS2905/coraalbers/ag/ag_data/gencode.v46.splice_sites_starts.feather'
HG38_SPLICE_END_PATH = '/users/PAS2905/coraalbers/ag/ag_data/gencode.v46.splice_sites_ends.feather'

HEART_UB = 'UBERON:0000948'
LV_UB = 'UBERON:0002084'
GENE = 'LMNA'

# # Flags to improve determinism.
# os.environ['XLA_FLAGS'] = ' '.join([
#     '--xla_gpu_deterministic_ops',
#     '--xla_gpu_enable_scatter_determinism_expander=True',
#     '--xla_gpu_enable_triton_gemm=False',
# ])
# # Increase GPU and CPU memory to reduce out of memory errors.
# os.environ['XLA_PYTHON_CLIENT_MEM_FRACTION'] = '0.9'

# gtf = pd.read_feather( 'https://storage.googleapis.com/alphagenome/reference/gencode/' 'hg38/gencode.v46.annotation.gtf.gz.feather' )

output_modalities = ['atac',	
    'cage',	
    'chip_histone',	
    'chip_tf',	
    'contact_maps',	
    'dnase',	
    'procap',	
    'rna_seq',	
    'splice_junctions',	
    'splice_site_usage',	
    'splice_sites']

requested_outputs = {dna_client.OutputType.ATAC,
        dna_client.OutputType.CAGE,
        dna_client.OutputType.DNASE,
        dna_client.OutputType.PROCAP,
        dna_client.OutputType.RNA_SEQ,
        dna_client.OutputType.SPLICE_SITES,
        dna_client.OutputType.SPLICE_SITE_USAGE,
        dna_client.OutputType.SPLICE_JUNCTIONS,
        dna_client.OutputType.CONTACT_MAPS,
        dna_client.OutputType.CHIP_HISTONE,
        dna_client.OutputType.CHIP_TF}

BED_FILE = sys.argv[1]
OUTPUT_FOLDER = sys.argv[2]


model = dna_model.create_from_huggingface( 
    'all_folds', 
    organism_settings={ 
        dna_model.Organism.HOMO_SAPIENS: dna_model.OrganismSettings( 
            fasta_path=HG38_FASTA_PATH, 
            gtf_feather_path=HG38_GTF_PATH, 
            splice_site_starts_feather_path=HG38_SPLICE_START_PATH, 
            splice_site_ends_feather_path=HG38_SPLICE_END_PATH, 
        ), dna_model.Organism.MUS_MUSCULUS: dna_model.OrganismSettings() } )

print('all_folds model initialized!')


def load_bed_intervals(bed_path, plot_interval=None):
    df = pd.read_csv(bed_path, sep='\t')
    intervals = [
        genome.Interval(
            chromosome=row.chrom,
            start=int(row.start),
            end=int(row.end),
            name=row['name'],
            info={'score': row.score},
        )
        for _, row in df.iterrows()
    ]
    num_intervals = len(intervals)
    if plot_interval is not None:
        intervals = [iv for iv in intervals if plot_interval.contains(iv)]
        print(f'loaded {num_intervals} intervals from bed file that fall within {plot_interval}!')
    else:
        print(f'loaded {num_intervals} intervals from bed file')
    return intervals, df


def tidy_scores_lv_chunked(variant_scores, ontology_curie, chunk_size=256):
    """Tidy ISM scores in chunks, keeping only tracks matching ontology_curie."""
    chunks = []
    n = len(variant_scores)
    for start in range(0, n, chunk_size):
        end = min(start + chunk_size, n)
        frames = []
        for i in range(start, end):
            for adata in variant_scores[i]:
                if adata is None or math.prod(adata.X.shape) == 0:
                    continue
                if 'ontology_curie' in adata.var.columns:
                    mask = adata.var['ontology_curie'] == ontology_curie
                    if not mask.any():
                        continue
                    adata = adata[:, mask]
                frame = variant_scorers.tidy_anndata(adata)
                if frame.empty:
                    continue
                if 'ontology_curie' in frame.columns:
                    frame = frame[frame['ontology_curie'] == ontology_curie]
                if not frame.empty:
                    frames.append(frame)
            variant_scores[i] = None  # free AnnDatas early
        if frames:
            chunks.append(pd.concat(frames, axis=0, ignore_index=True))
        del frames
        gc.collect()
    if not chunks:
        return pd.DataFrame()
    df = pd.concat(chunks, axis=0, ignore_index=True)
    last = ['raw_score'] + (['quantile_score'] if 'quantile_score' in df.columns else [])
    return df[[c for c in df.columns if c not in last] + last]


bed_intervals, bed_df = load_bed_intervals(BED_FILE)

scorers = dict(variant_scorers.RECOMMENDED_VARIANT_SCORERS)
scorers_drop = ['SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS', 'POLYADENYLATION']
scorers = {k: v for k, v in scorers.items() if k not in scorers_drop}
scorers = list(scorers.values())
print('scorers used')
print(*scorers, sep="\n")

os.makedirs(OUTPUT_FOLDER, exist_ok=True)

for interval_idx in range(len(bed_intervals)):
    
    ism_interval = bed_intervals[interval_idx]
    out_path = f"{OUTPUT_FOLDER}/var_scores_{ism_interval.name}.parquet"
    if os.path.exists(out_path):
        print(f'skipping {ism_interval.name}: already exists at {out_path}')
        continue

    print(f'interval {interval_idx + 1}/{len(bed_intervals)}: {ism_interval.name}')

    context = LMNA_INTERVAL.resize(dna_client.SEQUENCE_LENGTH_1MB)

    variant_scores = model.score_ism_variants( 
        interval=context, 
        ism_interval=ism_interval, 
        variant_scorers=scorers, 
        organism=dna_client.Organism.HOMO_SAPIENS, 
    )
    
    n_variants = len(variant_scores)
    n_scorers = len(variant_scores[0]) if n_variants else 0
    df = tidy_scores_lv_chunked(variant_scores, ontology_curie=LV_UB, chunk_size=256)
    del variant_scores
    gc.collect()

    # tidy_scores stores Variant/Interval objects; parquet needs strings.
    if 'variant_id' in df.columns:
        df['variant_id'] = df['variant_id'].map(str)
    if 'scored_interval' in df.columns:
        df['scored_interval'] = df['scored_interval'].map(str)
    
    print(f'number of variants analyzed: {n_variants}')
    print(f'number of scorers used per variant: {n_scorers}')
    print(f'number of rows in df: {len(df)}')

    df.to_parquet(out_path, compression="snappy")
    print(f'saved to {out_path}')

    del df
    gc.collect()










