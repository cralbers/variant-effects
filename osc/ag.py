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

if len(sys.argv) != 5:
    print("Usage: python ag.py <chrom> <target_start> <target_end> <region_num>")
    print("chrom formatted as chr#; target_start and target_end formatted as ###_###_### (use underscores every 3 digits)")
    sys.exit(1)

# initialize model
model = dna_model.create_from_huggingface( 
    'all_folds', 
    organism_settings={ 
        dna_model.Organism.HOMO_SAPIENS: dna_model.OrganismSettings( 
            fasta_path='hg38.fa', 
            gtf_feather_path='ag_data/gencode.v46.annotation.gtf.gz.feather', 
            splice_site_starts_feather_path='ag_data/gencode.v46.splice_sites_starts.feather', 
            splice_site_ends_feather_path='ag_data/gencode.v46.splice_sites_ends.feather', 
        ), dna_model.Organism.MUS_MUSCULUS: dna_model.OrganismSettings() } )



# DNA sequence to use as context when making predictions.
chrom = sys.argv[1]
target_start = int(sys.argv[2])
target_end = int(sys.argv[3])
region_num = sys.argv[4]

# define intervals
ism_interval = genome.Interval(chrom, target_start, target_end) 
# sequence_interval = ism_interval.resize(32768) 
lmna_context = ism_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

sj_variant_scorer = variant_scorers.RECOMMENDED_VARIANT_SCORERS['SPLICE_JUNCTIONS']

variant_scores = model.score_ism_variants( 
    interval=lmna_context, 
    ism_interval=ism_interval, 
    variant_scorers=[sj_variant_scorer], 
    organism=dna_client.Organism.HOMO_SAPIENS, 
)

variant_scores.to_pickle(f"variant_scores_region{region_num}.pkl")

def extract_splice_scores_fixed(adata): 
    var = adata.var 
    mask = (var['ontology_curie'] == 'UBERON:0000948') & (var['Assay title'] == 'total RNA-seq') 
    if adata.X.shape[0] == 0: 
        return float('nan') 
    values = adata.X[0, mask] 
    if values.size == 0: 
        return float('nan') 
    return float(values.max())
    
rows = [] 
for vs in variant_scores: 
    variant = vs[0].uns['variant'] 
    score = extract_splice_scores_fixed(vs[0]) 
    rows.append({ 
        'position': variant.position, 
        'ref': variant.reference_bases, 
        'alt': variant.alternate_bases, 
        'delta_score': score, 
    })



ism_df = pd.DataFrame(rows) 
ism_df.to_csv(f'lmna_ism_out/lmna_ism_region{region_num}.csv', index=False) 
print(f"Total rows: {len(ism_df)}") 
print(f"Non-NaN: {ism_df['delta_score'].notna().sum()}") 
print(f"NaN: {ism_df['delta_score'].isna().sum()}") 
print("Done!")







