#!/bin/bash
#SBATCH -J fst_climate_15
#SBATCH -p barbun
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/fst_results/fst_climate_15_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/fst_results/fst_climate_15_%j.err


cd /arf/scratch/ssenkal/selin/sekanslar/all_years_26/fst_results

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

Rscript PCA_ancova_iklim_fst_15var.R

echo ">>> Tamamlandi: $(date)"
