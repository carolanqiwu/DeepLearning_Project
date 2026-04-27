#!/bin/bash
#SBATCH --job-name="GBM_hic_extract"
#SBATCH --partition=lesliec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=0-00:30:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wuc10@mskcc.org
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

#RUN this to sbatch slurm array
#N=$(find /data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/hichip \
#    -maxdepth 1 -type d -name "s*" | wc -l)

#sbatch --array=0-$(($N-1))%2 extract_hic_txt.sh

mkdir -p logs
source /data1/lesliec/carolw/my_envs/archR/etc/profile.d/conda.sh
conda activate /data1/lesliec/carolw/my_envs/chromafold

JAR_DIR=/data1/lesliec/carolw/scripts/hicdc/juicer_tools.1.9.9_jcuda.0.8.jar
hichip_dir=/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan/hic_GSC

# Get list of GBMx-* folders
mapfile -t Folder_dir < <(find "${hichip_dir}" -maxdepth 1 -type d -name "s*" | sort)

# Select folder for this array task
FOLDER="${Folder_dir[$SLURM_ARRAY_TASK_ID]}"
BASENAME=$(basename "$FOLDER")

PREFIX="hg38_10kb_GATC_GANTC_FDR_05"

zvalue_hic=${FOLDER}/${BASENAME}_${PREFIX}_zvalue.hic
qvalue_hic=${FOLDER}/${BASENAME}_${PREFIX}_qvalue.hic
out_zvalue=${FOLDER}/zvalue
out_qvalue=${FOLDER}/qvalue

mkdir -p $out_zvalue
mkdir -p $out_qvalue

for i in {1..22}
do

java -jar  $JAR_DIR dump observed NONE $zvalue_hic $i $i BP 10000 $out_zvalue/chr"$i"_raw.txt

java -jar  $JAR_DIR dump observed NONE $qvalue_hic $i $i BP 10000 $out_qvalue/chr"$i"_raw.txt
echo $i

done
