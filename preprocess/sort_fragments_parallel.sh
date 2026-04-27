#!/bin/bash
#SBATCH --job-name="GSC_sort_fragments"
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

BASE_DIR="/data1/lesliec/carolw/data_raw/GSC"
SORTED_DIR="/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/fragments_sorted"

# find all fragments.tsv.gz files under GSC
find "${BASE_DIR}" -type f -name "fragments.tsv.gz" | sort | while read -r FRAG_FILE; do

    # get subfolder name (parent directory)
    SUBFOLDER=$(basename "$(dirname "${FRAG_FILE}")")

    # set data prefix as GBMx_<subfolder name>
    DATA_PREFIX="GBMx_${SUBFOLDER}"

    echo "Processing ${DATA_PREFIX}"

    # Sorted output path
    SORTED_TSV="${SORTED_DIR}/${DATA_PREFIX}_bgz_sorted.tsv"

    # sort + bgzip
    gunzip -c "${FRAG_FILE}" \
        | sortBed -i stdin \
        > "${SORTED_TSV}"

    htsfile "${SORTED_TSV}"
    bgzip -f "${SORTED_TSV}"

done
