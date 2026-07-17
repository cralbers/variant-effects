#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --partition=nextgen
#SBATCH --job-name=lmna_reg1
#SBATCH --mem=64gb
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1


module load cuda/12.9.1
module load miniconda3/24.1.2-py310


source activate py311

nvidia-smi
module list
echo $CUDA_VISIBLE_DEVICES
echo $LD_LIBRARY_PATH
python -c "import jax; print(jax.devices())"


python ag.py chr1 156104711 156114711 1