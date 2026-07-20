#!/usr/bin/env nextflow


params.vat_path = "gs://vwb-aou-datasets-controlled/v9/wgs/short_read/snpindel/aux/vat/*.gz"
params.gene_symbol = "LMNA"
params.gene_id = "ENSG00000160789"   
params.output_ht = "${WORKSPACE_BUCKET}/hail_files/lmna_vat.ht"
params.google_project = "${GOOGLE_CLOUD_PROJECT}"
params.docker = "${ARTIFACT_REGISTRY_DOCKER_REPO}/hailgenetics/hail:0.2.128"
params.cpus = 4
params.memory = "20GB"
params.disk = "200.GB"

process FilterLmnaVat {
    container params.docker
    cpus params.cpus
    memory params.memory
    disk params.disk

    // Publish the log; the .ht lives on GCS (written by Hail), not as a Nextflow path output
    publishDir "results", mode: "copy"

    input:
    path python_script
    val vat_path
    val output_ht
    val gene_id
    val gene_symbol
    val google_project

    output:
    path "run.log"
    val output_ht, emit: hail_table

    script:
    """
    set -euo pipefail
    export GOOGLE_CLOUD_PROJECT="${google_project}"
    run Setting_Env_Variables_p2.py


    python3 ${python_script} \\
      --vat-path "${vat_path}" \\
      --output-ht "${output_ht}" \\
      --gene-id "${gene_id}" \\
      --gene-symbol "${gene_symbol}" \\
      --google-project "${google_project}" \\
      --overwrite \\
      2>&1 | tee run.log
    """
}

workflow {
    if( !params.output_ht || !params.google_project ) {
        error "Set --output_ht and --google_project (gs://... and GCP project id)"
    }

    FilterLmnaVat(
        file("${projectDir}/filter_lmna_vat.py"),
        params.vat_path,
        params.output_ht,
        params.gene_id,
        params.gene_symbol,
        params.google_project
    )
}

workflow.onComplete {
    println( workflow.success
        ? "OK: wrote ${params.output_ht}"
        : "Failed: ${workflow.errorReport}"
    )
}
