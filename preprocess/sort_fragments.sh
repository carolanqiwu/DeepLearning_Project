#!/bin/bash
#SBATCH --job-name="GBM_Proneural_sort_fragments"
#SBATCH --partition=lesliec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wuc10@mskcc.org
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

mkdir -p logs
source /admin/software/miniforge3/etc/profile.d/conda.sh
conda activate /data1/lesliec/carolw/my_envs/chromafold

SORTED_DIR="/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/archr_data_subsampled_nonimmune"
FRAG_FILE="${SORTED_DIR}/NonImmune_Proneural_2k_fragments.tsv.gz"
DATA_PREFIX="NonImmune_Proneural_2k"

# Sorted output path
SORTED_TSV="${SORTED_DIR}/${DATA_PREFIX}_bgz_sorted.tsv"
# sort + bgzip
gunzip -c "${FRAG_FILE}" \
	| sort -k 1,1 -k2,2n \
        > "${SORTED_TSV}"

htsfile "${SORTED_TSV}"
bgzip -f "${SORTED_TSV}"

