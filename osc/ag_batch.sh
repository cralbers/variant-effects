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
#SBATCH --gpus-per-task=1

# load the necessary modules
module load cuda/12.9.1
module load miniconda3/24.1.2-py310

# activate the alphagenome conda environment
source activate py311

# these lines are used to check the environment (modules loaded)
# and to print the GPU information (JAX devices)
# to confirm that the GPU is being detected and used
nvidia-smi
module list
echo $CUDA_VISIBLE_DEVICES
echo $LD_LIBRARY_PATH
python -c "import jax; print(jax.devices())"


## run script below