#!/bin/bash
#SBATCH --job-name="GBM_CLS_chromafold_predict"
#SBATCH --partition=lesliec
#SBATCH --gres=gpu:a100:1
#SBATCH --mem=150G
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wuc10@mskcc.org
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/chromafold_%A_%a.out
#SBATCH --error=logs/chromafold_%A_%a.err

#RUN this as slurm array
#N=$(ls /data1/lesliec/carolw/projects/chromafold/preprocess/BRCA-chan/archr_data | wc -l)
#sbatch --array=0-$(($N-1))%1 run_inference.slurm
mkdir -p logs

#module load R/4.3.0
source /admin/software/miniforge3/etc/profile.d/conda.sh
conda activate /data1/lesliec/carolw/my_envs/chromafold

set -euo pipefail

SCRIPT_DIR="/data1/lesliec/carolw/repos/ChromaFold/chromafold"
DATA_PATH="/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan"
SAVE_PATH="/data1/lesliec/carolw/projects/chromafold/predictions/GBM-chan/GBMx_NonImmune_Classical_2k"
MODEL_PATH="/data1/lesliec/carolw/repos/ChromaFold/checkpoints/chromafold_CTCFmotif.pth.tar"
GENOME_ASSEMBLY="hg38"
#TEST_CHROM=( {1..22} X )
TEST_CHROM=({1,4,7,19})
PREFIX="GBMx_NonImmune_Classical_2k"
mkdir -p "${SAVE_PATH}"
# ===============================
# Collect ArchR subfolders
# ===============================

python inference.py --data-path "${DATA_PATH}" \
--model-path "${MODEL_PATH}" \
--save-path "${SAVE_PATH}" \
-ct "${PREFIX}" \
--test-chrom "${TEST_CHROM[@]}" \
-offset -2000000 \
--genome "${GENOME_ASSEMBLY}"


