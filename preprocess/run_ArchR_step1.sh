#!/bin/bash
#SBATCH --job-name="GBM_ArchR_prep1"
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
set -eo pipefail 

source /admin/software/miniforge3/etc/profile.d/conda.sh
conda activate /data1/lesliec/carolw/my_envs/archR

set -euo pipefail
#RUN this to sbatch script:
#N=$(ls /data1/lesliec/carolw/data_raw/scATAC-chan/Cancer_scATACseq_data/scATAC_BRCA_*.fragments.tsv.gz | wc -l)
#sbatch --array=0-$(($N-1)) run_ArchR.sh

SAVE_LOC="/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan"
FRAG_LOC="/data1/lesliec/carolw/data_raw/scATAC-chan/Cancer_scATACseq_data"
GENOME_ASSEMBLY="hg38"
SORTED_DIR="${SAVE_LOC}/fragments_sorted"
mkdir -p "${SORTED_DIR}"

FILES=(${FRAG_LOC}/scATAC_GBMx_*.fragments.tsv.gz)
FRAG_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"
if [[ -z "${FRAG_FILE}" ]]; then
  echo "No fragment file for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Fragment file: ${FRAG_FILE}"

# extract prefix: BRCA_XX
DATA_PREFIX=$(basename "$FRAG_FILE" | sed -E 's/^scATAC_(GBMx_[^.]*)\.fragments\.tsv\.gz/\1/')
FRAG_FILE_PREFIX="${DATA_PREFIX}"

# ArchR output location
ARCHR_LOC="${SAVE_LOC}/archr_data/${DATA_PREFIX}"
mkdir -p "${ARCHR_LOC}"

# use linux sort
SORTED_TSV="${SORTED_DIR}/${DATA_PREFIX}_bgz_sorted.tsv"

gunzip -c "${FRAG_FILE}" | sort -k1,1 -k2,2n > "${SORTED_TSV}"
htsfile "${SORTED_TSV}"
bgzip -f "${SORTED_TSV}"

#rm "${FRAG_LOC}"/"${FRAG_FILE_PREFIX}_bgz_sorted.tsv.gz.tbi" # remove previously calculated .tbi file

# create folders to store ChromaFold input
mkdir -p "${SAVE_LOC}/atac" "${SAVE_LOC}/dna" "${SAVE_LOC}/predictions"
# copy CTCF motif data for each genome assembly
cp -n /data1/lesliec/carolw/repos/ChromaFold/hg38_ctcf_motif_score.p "${SAVE_LOC}"/dna/

# run R to create LSI file using ArchR
Rscript ArchR_preparation.R \
    "${DATA_PREFIX}" \
    "${ARCHR_LOC}" \
    "${SORTED_TSV}.gz" \
    "${GENOME_ASSEMBLY}" \
    1365
   
echo "Finished ${DATA_PREFIX}"

