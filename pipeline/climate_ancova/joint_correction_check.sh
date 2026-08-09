#!/bin/bash
#SBATCH -p barbun
#SBATCH -A ssenkal
#SBATCH -J joint_correction_check
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=32G
#SBATCH -t 01:00:00
#SBATCH --qos=normal
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/joint_corr_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/joint_corr_%j.err

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

Rscript /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/joint_correction_check.R
