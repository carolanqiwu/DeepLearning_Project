#!/bin/bash
#SBATCH --job-name="chromafold_vis"
#SBATCH --partition=lesliec
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=128G
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=wuc10@mskcc.org


# Load conda
source /admin/software/anaconda/1.11.1/etc/profile.d/conda.sh
/admin/software/anaconda/1.11.1/envs/jupyter/bin/python vis_eval.py


