#!/bin/bash
#SBATCH -p hamsi
#SBATCH -A ssenkal
#SBATCH -J inversion_freq
#SBATCH --nodes=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --qos=normal
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/INVERSION_ANALYSIS/inv_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/INVERSION_ANALYSIS/inv_%j.err

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

python3 /arf/scratch/ssenkal/selin/sekanslar/all_years_26/INVERSION_ANALYSIS/extract_inversion_freqs.py
Rscript /arf/scratch/ssenkal/selin/sekanslar/all_years_26/INVERSION_ANALYSIS/inversion_summary_plot.R
