#!/bin/bash
#SBATCH -J fst_corrected
#SBATCH -p barbun
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH -t 24:00:00
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/fst_results/fst_corrected_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/fst_results/fst_corrected_%j.err

PROGRAMS="/arf/scratch/ssenkal/selin/programs/popoolation2_1201"
OUTDIR="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/fst_results"
SYNC="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/mega_all_years.sync"

mkdir -p $OUTDIR
source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate align_env

echo ">>> Popoolation2 Sliding Window Fst (corrected pool-size) Başlıyor..."

# Sıra (BAM_LIST ile doğrulandı): 2015_06,08,09,10, 2016_06,08,09,10, 2017_06,10,
#       2018_06,08,09,10, 2019_06,08,09,10, 2020_06,08,09,10, 2021_06,08,09,10
# 17. pozisyon = 2019_09 = 30 birey -> pool-size 60
# Diğer tüm 25 kohort = 40 birey -> pool-size 80
POOLSIZES="80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,80,60,80,80,80,80,80,80,80,80,80"

perl $PROGRAMS/fst-sliding.pl \
    --input $SYNC \
    --output $OUTDIR/mega_all_years_26_sliding_fst_corrected.txt \
    --pool-size $POOLSIZES \
    --min-count 2 \
    --min-coverage 15 \
    --max-coverage 400 \
    --window-size 5000 \
    --step-size 5000 \
    --suppress-noninformative

echo ">>> İşlem Tamamlandı: $(date)"
