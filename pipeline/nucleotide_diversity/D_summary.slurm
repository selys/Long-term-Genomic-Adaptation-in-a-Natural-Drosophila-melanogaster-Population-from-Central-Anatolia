#!/bin/bash
#SBATCH -J D_summary_ps80
#SBATCH -p barbun
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task=40
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=/arf/scratch/ssenkal/selin/sekanslar/all_years_26/PI_ANALYSIS_1k_ps80/verify_ps80/D_summary_%j.out
#SBATCH --error=/arf/scratch/ssenkal/selin/sekanslar/all_years_26/PI_ANALYSIS_1k_ps80/verify_ps80/D_summary_%j.err

cd /arf/scratch/ssenkal/selin/sekanslar/all_years_26/PI_ANALYSIS_1k_ps80/verify_ps80

echo -e "sample\tmean_D" > D_summary_ps80_corrected.tsv

for f in *_D_ps80.txt; do
    sample=$(basename "$f" _D_ps80.txt)
    mean_d=$(awk '{sum+=$5; count++} END {if(count>0) print sum/count; else print "NA"}' "$f")
    echo -e "${sample}\t${mean_d}" >> D_summary_ps80_corrected.tsv
done

echo ">>> D_summary_ps80_corrected.tsv olusturuldu: $(date)"
cat D_summary_ps80_corrected.tsv
