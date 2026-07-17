#!/bin/bash
#SBATCH --account=pas2905
#SBATCH --partition=nextgen
#SBATCH --job-name=ism_reg_regions
#SBATCH --time=40:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --output=R-%x.%j.out
#SBATCH --error=R-%x.%j.err
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=10
#SBATCH --gpus-per-task=1   # or 2 if you really need two local GPUs


module load cuda/12.9.1
module load miniconda3/24.1.2-py310


source activate py311

nvidia-smi
module list
echo $CUDA_VISIBLE_DEVICES
echo $LD_LIBRARY_PATH
python -c "import jax; print(jax.devices())"

export XLA_PYTHON_CLIENT_MEM_FRACTION=0.9
export XLA_FLAGS='--xla_gpu_deterministic_ops --xla_gpu_enable_scatter_determinism_expander=True --xla_gpu_enable_triton_gemm=False'

python run_ism_regulatory_regions.py 'outputs/predicted_reg_regions_post_combinedtracks.bed' 'outputs/post_var_scores'