#!/bin/bash
#SBATCH --job-name="GBM_Classical_ArchR_prep2"
#SBATCH --partition=lesliec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=150G
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wuc10@mskcc.org


source /admin/software/miniforge3/etc/profile.d/conda.sh
conda activate /data1/lesliec/carolw/my_envs/chromafold

SCRIPT_DIR="/data1/lesliec/carolw/repos/ChromaFold/preprocessing_pipeline"
SAVE_LOC="/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan"
GENOME_ASSEMBLY="hg38"

ARCHR_LOC="${SAVE_LOC}/archr_data_subsampled_nonimmune"
SORTED_FRAG="${ARCHR_LOC}/NonImmune_Classical_2k_bgz_sorted.tsv.gz"

(echo '"","x"'; tail -n +2 "${ARCHR_LOC}"/Classical_barcodes_2k.csv | nl -w1 -s',' -v1 | sed 's/\t/,/') > "${ARCHR_LOC}"/Classical_2k_numbered_barcodes.csv
# step 2 calculate tile files using ArchR
python "${SCRIPT_DIR}"/scATAC_preparation.py \
--cell_type_prefix "GBMx_NonImmune_Classical_2k" \
--fragment_file  "${SORTED_FRAG}" \
--barcode_file "${ARCHR_LOC}"/Classical_2k_numbered_barcodes.csv \
--lsi_file "${ARCHR_LOC}"/Classical_lsi_2k.csv \
--genome_assembly "${GENOME_ASSEMBLY}" \
--save_path "${SAVE_LOC}"

echo "Finished scATAC step 2"
