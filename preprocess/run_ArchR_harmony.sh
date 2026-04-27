#!/bin/bash
#SBATCH --job-name="GBM_ArchR_prep1"
#SBATCH --partition=lesliec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wuc10@mskcc.org
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

mkdir -p logs

source /admin/software/miniforge3/etc/profile.d/conda.sh
conda activate /data1/lesliec/carolw/my_envs/archr_2

set -euo pipefail

SAVE_LOC="/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan"
mkdir -p "${SAVE_LOC}"

FRAG_LOC="${SAVE_LOC}/fragments_sorted"
GENOME_ASSEMBLY="hg38"
SUBSAMPLE=FALSE
COHORT="GBMx"
#CELL_NUMBER=1000
# ArchR output location
ARCHR_LOC="${SAVE_LOC}/archr_data_allcells"
mkdir -p "${ARCHR_LOC}"

# create folders to store ChromaFold input
mkdir -p "${SAVE_LOC}/atac" "${SAVE_LOC}/dna" 
# copy CTCF motif data for each genome assembly
cp -n /data1/lesliec/carolw/repos/ChromaFold/hg38_ctcf_motif_score.p "${SAVE_LOC}"/dna/

# run R to create LSI file using ArchR
Rscript ArchR_preparation_harmony.R \
    "${ARCHR_LOC}" \
    "${FRAG_LOC}" \
    "${GENOME_ASSEMBLY}" \
    "${SUBSAMPLE}" \
    "${COHORT}" 

echo "Finished ArchR preparation with harmony"

