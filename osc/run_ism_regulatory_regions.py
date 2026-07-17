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

if len(sys.argv) != 3:
    print("Usage: python run_ism_regulatory_regions.py <bed_file_path> <output_folder>")
    print("")
    sys.exit(1)

LMNA_START = 156_114_711
LMNA_END = 156_140_081
gene_symbol = "LMNA"
LMNA_INTERVAL = genome.Interval('chr1', 156_114_711, 156_140_081)


BASE_PATH = '/users/PAS2905/coraalbers/'
AG_DATA_PATH = '/users/PAS2905/coraalbers/ag/ag_data/'

HG38_FASTA_PATH = '/users/PAS2905/coraalbers/ag/hg38.fa'
HG38_GTF_PATH = '/users/PAS2905/coraalbers/ag/ag_data/gencode.v46.annotation.gtf.gz.feather'
HG38_SPLICE_START_PATH = '/users/PAS2905/coraalbers/ag/ag_data/gencode.v46.splice_sites_starts.feather'
HG38_SPLICE_END_PATH = '/users/PAS2905/coraalbers/ag/ag_data/gencode.v46.splice_sites_ends.feather'

HEART_UB = 'UBERON:0000948'
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


bed_intervals, bed_df = load_bed_intervals(BED_FILE)

scorers = dict(variant_scorers.RECOMMENDED_VARIANT_SCORERS)
scorers_drop = ['SPLICE_SITES', 'SPLICE_SITE_USAGE', 'SPLICE_JUNCTIONS', 'POLYADENYLATION']
scorers = {k: v for k, v in scorers.items() if k not in scorers_drop}
scorers = list(scorers.values())
print('scorers used')
print(*scorers, sep="\n")
    

for i in range(len(bed_intervals)):
    
    ism_interval = bed_intervals[i]
    print(f'interval name: {ism_interval.name}')

    context = ism_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

    variant_scores = model.score_ism_variants( 
        interval=context, 
        ism_interval=ism_interval, 
        variant_scorers=scorers, 
        organism=dna_client.Organism.HOMO_SAPIENS, 
    )

    print(f'number of variants analyzed: {len(variant_scores)}')
    print(f'number of scorers used per variant: {len(variant_scores[0])}')
    with open(f"{OUTPUT_FOLDER}/var_scores_{ism_interval.name}.pkl", "wb") as f:
        pickle.dump(variant_scores, f)










