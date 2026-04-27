#!/bin/bash
#SBATCH --job-name="GBMx_9976F952_ArchR_cluster"
#SBATCH --partition=lesliec
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wuc10@mskcc.org
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

mkdir -p logs

source /admin/software/miniforge3/etc/profile.d/conda.sh
conda activate /data1/lesliec/carolw/my_envs/archr_2

base_path="/data1/lesliec/carolw/projects/chromafold/preprocess/GBM-chan"
archr_path="${base_path}/1_3k/archr_data/GBMx_9976F952_23A5_431A_A431_01E544324A26_X010_S05_B1_T1"
out_path="${base_path}/GBMx_9976F952"
genome_assembly="hg38"

mkdir -p "${out_path}"
Rscript ArchR_cluster_analysis.R ${archr_path} ${out_path} ${genome_assembly}
