#!/bin/bash
#SBATCH -p barbun
#SBATCH -A ssenkal
#SBATCH -J ANCOVA_no2015_09
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH -t 24:00:00
#SBATCH --qos=normal
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/ancova_no2015_09_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/ancova_no2015_09_%j.err

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

echo ">>> Job started: $(date)"
Rscript /arf/scratch/ssenkal/selin/sekanslar/snp_climate_ancova_no2015_09.R
echo ">>> Job finished: $(date)"
