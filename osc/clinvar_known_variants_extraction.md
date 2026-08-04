# Get noncoding variants

GTF=/users/PAS2905/coraalbers/ag/ag_data/gencode.v46.annotation.gtf
GENOME=hg38.chrom.sizes   # chrom\tlength

1. download clinvar dataset
```bash
curl -O https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260728.vcf.gz
```

2. Get CDS intervals (GTF is 1-based inclusive → BED is 0-based half-open)

```bash
awk 'BEGIN{OFS="\t"} $3=="CDS" {print $1,$4-1,$5}' /users/PAS2905/coraalbers/ag/ag_data/gencode.v46.annotation.gtf | bedtools sort -i - | bedtools merge -i - > /users/PAS2905/coraalbers/ag/ag_data/cds.bed
```

**if needed, get chromosome sizes from fasta index file (only needs to be done once)**

```bash
cut -f1,2 /users/PAS2905/coraalbers/ag/hg38.fa.fai > hg38.chrom.sizes
```

3. Extract noncoding regions (genome complement of CDS)

**needs to be run in folder with fasta index file and chromosome sizes file**

```bash
bedtools complement -i /users/PAS2905/coraalbers/ag/ag_data/cds.bed -g /users/PAS2905/coraalbers/ag/hg38.chrom.sizes > /users/PAS2905/coraalbers/ag/ag_data/noncoding.bed
```

4. remove chr annotation to match vcf file chromosome naming

```bash
sed 's/^chr//' noncoding.bed > noncoding_num_chr.bed
```

5. run clinvar variant extraction notebook to generate vcfs of variants within 500kb region of gene
6. use bedtools intersect to get variants that overlap with noncoding regions

```bash
bedtools intersect -wo -a /users/PAS2905/coraalbers/ag/variant-effects/osc/outputs/clinvar_LMNA.PLP.vcf -b /users/PAS2905/coraalbers/ag/ag_data/noncoding_num_chr.bed > outputs/plp_with_nc.bed
bedtools intersect -wo -a /users/PAS2905/coraalbers/ag/variant-effects/osc/outputs/clinvar_LMNA.BLB.vcf -b /users/PAS2905/coraalbers/ag/ag_data/noncoding_num_chr.bed > outputs/blb_with_nc.bed
```

