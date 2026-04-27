#!/bin/bash
#SBATCH --job-name="BRCA_merged_ArchR_cluster"
#SBATCH --partition=lesliec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=64G
#SBATCH --time=1-00:00:00
#SBATCH --array=0-9%10
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wuc10@mskcc.org
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

mkdir -p logs

source /admin/software/miniforge3/etc/profile.d/conda.sh
conda activate /data1/lesliec/carolw/my_envs/archr_2

base_path="/data1/lesliec/carolw/projects/chromafold/preprocess/BRCA-chan"
genome_assembly="hg38"

# Build array of all subfolders
archr_paths=($(ls -d "${base_path}/1_3k/archr_data"/*/))

# Get the folder for this array task
archr_path="${archr_paths[$SLURM_ARRAY_TASK_ID]}"
subfolder=$(basename "${archr_path}")
prefix=$(echo "${subfolder}" | cut -d'_' -f1,2)
out_path="${base_path}/1_3k/${prefix}_UMAP"

mkdir -p "${out_path}"
mkdir -p logs

Rscript ArchR_cluster_analysis.R ${archr_path} ${out_path} ${genome_assembly}
