#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --partition=nextgen
#SBATCH --job-name=lmna_ism_1mb
#SBATCH --mem=128gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=logs/R-%x.%j.out
#SBATCH --error=logs/R-%x.%j.err


cd /users/PAS2905/coraalbers/ag/variant-effects/osc

# load the necessary modules
ml miniconda3/24.1.2-py310 cuda/12.9.1

# activate the alphagenome conda environment
conda activate py311

python api_based_ism_recommended_scorers.py
     
