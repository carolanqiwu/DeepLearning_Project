#!/bin/bash
#SBATCH --job-name="merge_fragments"
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
conda activate /data1/lesliec/carolw/my_envs/archr_2

set -euo pipefail

IN_DIR="/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/fragments_sorted"
OUT_DIR="/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/fragments_sorted/merged"

mkdir -p "${OUT_DIR}"

for f in "${IN_DIR}"/*.tsv.gz; do
  echo "Processing $f"

  # Extract prefix from filename
  # Keep only BRCA_08499A64 (first two underscore-separated fields)
  prefix="GBMx_$(basename "$f" | cut -d'_' -f2)"

  out="${OUT_DIR}/${prefix}_fragments.tsv"
  
  gzip -dc "$f" | \
    awk -v pfx="${prefix}#" 'BEGIN{FS=OFS="\t"} {print $1,$2,$3,pfx$4,$5}'  \
    > "$out"

done

echo "Sorting per-sample files"
for f in "${OUT_DIR}"/GBMx_*_fragments.tsv; do
  sort -T "${OUT_DIR}" --buffer-size=2G -k1,1V -k2,2n "$f" -o "$f"
done


echo "Merging fragment files..."
sort -T "${OUT_DIR}" -k1,1V -k2,2n \
  "${OUT_DIR}"/GBMx_*_fragments.tsv \
  > "${OUT_DIR}/GBMx_merged_fragments.tsv"

#cat "${TMP_FILES[@]}" | sort -T "${OUT_DIR}" -k1,1V -k2,2n > "${OUT_DIR}/GBMx_merged_fragments.tsv"
#cat "${OUT_DIR}"/BRCA_*_fragments.tsv | sort -T "${OUT_DIR}" -k1,1V -k2,2n > "${OUT_DIR}/BRCA_merged_fragments.tsv"

echo "Compressing merged fragments..."
bgzip -@ 4 "${OUT_DIR}/GBMx_merged_fragments.tsv"

echo "Indexing with tabix..."
tabix -p bed "${OUT_DIR}/GBMx_merged_fragments.tsv.gz"

#echo "Cleaning up intermediate files..."
#rm "${TMP_FILES[@]}"


echo "Finished merging"

