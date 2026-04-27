#!/bin/bash
#SBATCH --job-name="BRCA_ArchR_prep2"
#SBATCH --partition=lesliec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wuc10@mskcc.org

#RUN this to sbatch script
#N=$(ls /data1/lesliec/carolw/data_raw/scATAC-chan/Cancer_scATACseq_data/scATAC_BRCA_*.fragments.tsv.gz | wc -l)
#sbatch --array=0-$(($N-1)) run_ArchR_step2.sh

source /admin/software/miniforge3/etc/profile.d/conda.sh
conda activate /data1/lesliec/carolw/my_envs/chromafold

SCRIPT_DIR="/data1/lesliec/carolw/repos/ChromaFold/preprocessing_pipeline"
SAVE_LOC="/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA-chan"
FRAG_LOC="/data1/lesliec/carolw/data_raw/scATAC-chan/Cancer_scATACseq_data"
GENOME_ASSEMBLY="hg38"

# discover fragment files (raw names define sample set)
FILES=(${FRAG_LOC}/scATAC_BRCA_*.fragments.tsv.gz)
FRAG_FILE="${FILES[$SLURM_ARRAY_TASK_ID]}"

if [[ -z "${FRAG_FILE}" ]]; then
  echo "No fragment file for task ${SLURM_ARRAY_TASK_ID}"
  exit 1
fi

echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "Fragment file: ${FRAG_FILE}"

# extract BLCA_XX
DATA_PREFIX=$(basename "${FRAG_FILE}" | sed -E 's/^scATAC_(BRCA_[^.]*)\.fragments\.tsv\.gz/\1/')

ARCHR_LOC="${SAVE_LOC}/archr_data/${DATA_PREFIX}"
SORTED_DIR="${SAVE_LOC}/fragments_sorted"
SORTED_FRAG="${SORTED_DIR}/${DATA_PREFIX}_bgz_sorted.tsv.gz"

# sanity checks
if [[ ! -f "${SORTED_FRAG}" ]]; then
  echo "ERROR: sorted fragments not found: ${SORTED_FRAG}"
  exit 1
fi

if [[ ! -f "${ARCHR_LOC}/archr_filtered_barcode.csv" ]]; then
  echo "ERROR: missing ArchR barcode file for ${DATA_PREFIX}"
  exit 1
fi

if [[ ! -f "${ARCHR_LOC}/archr_filtered_lsi.csv" ]]; then
  echo "ERROR: missing ArchR LSI file for ${DATA_PREFIX}"
  exit 1
fi

# step 2 calculate tile files using ArchR
python "${SCRIPT_DIR}"/scATAC_preparation.py \
--cell_type_prefix "${DATA_PREFIX}" \
--fragment_file  "${SORTED_FRAG}" \
--barcode_file "${ARCHR_LOC}"/archr_filtered_barcode.csv \
--lsi_file "${ARCHR_LOC}"/archr_filtered_lsi.csv \
--genome_assembly "${GENOME_ASSEMBLY}" \
--save_path "${SAVE_LOC}"

echo "Finished scATAC step 2 for ${DATA_PREFIX}"
