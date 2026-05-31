====================================================================
1. BIOINFORMATICS PIPELINE (TRUBA BASH SCRIPTS)
====================================================================
# This section contains the Pool-Seq quality control, alignment, 
# sorting, and variant calling pipelines executed on the TRUBA cluster.


#!/bin/bash
#SBATCH -J bwa_align
#SBATCH -p barbun
#SBATCH --array=1-40%4
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH -o logs/bwa_%A_%a.out
#SBATCH -e logs/bwa_%A_%a.err

set -euo pipefail

#####################
# CONDA
#####################
source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate align_env

#####################
# PATHS
#####################
BASE=/arf/scratch/ssenkal/selin/sekanslar
TRIM_DIR=$BASE/qc_fastp_all
BAM_DIR=$BASE/bam
REF=$BASE/RAW_DIR/dmel_reference_genome/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna

mkdir -p $BAM_DIR logs

#####################
# GET SRR
#####################
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $TRIM_DIR/srr_list.txt)

R1=$TRIM_DIR/$SRR/${SRR}_1.trim.fastq.gz
R2=$TRIM_DIR/$SRR/${SRR}_2.trim.fastq.gz

#####################
# ALIGNMENT
#####################
bwa mem -t 8 $REF $R1 $R2 \
| samtools view -bh -q 20 -F 0x100 \
| samtools sort -@ 4 -o $BAM_DIR/${SRR}.sorted.bam

samtools index $BAM_DIR/${SRR}.sorted.bam


#!/bin/bash
#SBATCH -p orfoz
#SBATCH -A ssenkal
#SBATCH -J BAM_QC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=32G
#SBATCH -t 04:00:00
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/LOGS/qc_%j_%x.out

source $(conda info --base)/etc/profile.d/conda.sh
conda activate drosophila_env

SAMPLE=$1
OUT_DIR="/arf/scratch/ssenkal/selin/sekanslar/MAPPED_DATA"
QC_DIR="/arf/scratch/ssenkal/selin/sekanslar/QC_REPORTS"

mkdir -p $QC_DIR

echo "--- QC Basliyor: $SAMPLE ---"

# 1. Flagstat (Eşleşme oranları)
samtools flagstat ${OUT_DIR}/${SAMPLE}_final.bam > ${QC_DIR}/${SAMPLE}_flagstat.txt

# 2. Coverage ve Ortalama Derinlik (Hızlı özet)
samtools coverage ${OUT_DIR}/${SAMPLE}_final.bam > ${QC_DIR}/${SAMPLE}_coverage.txt

# 3. Depth (Daha detaylı analiz için her bazın derinliği)
# Bu komut ağır olabilir, sadece genel ortalamayı almak için 'stats' kullanalım
samtools stats ${OUT_DIR}/${SAMPLE}_final.bam | grep ^SN > ${QC_DIR}/${SAMPLE}_stats.txt

echo "--- $SAMPLE QC Tamamlandi ---"

————MPILEUP and ANNOTATION STEPS——
#SBATCH -p hamsi
#SBATCH -A ssenkal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH -t 36:00:00
#SBATCH -J mega_full_genom_26
#SBATCH -o mega_sync_%j.out
#SBATCH -e mega_sync_%j.err
#SBATCH --nodelist=hamsi11

# Ortam ve Yollar
source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

YEAR="mega_all_years"
BASE_DIR="/arf/scratch/ssenkal/selin/sekanslar"
REF="$BASE_DIR/RAW_DIR/dmel_reference_genome/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna"
BAM_LIST="$BASE_DIR/all_years_jun_aug_sep_oct_26.txt"
OUT_DIR="$BASE_DIR/all_years_26"
mkdir -p $OUT_DIR

# BAM dosyalarını oku
BAM_FILES=$(cat $BAM_LIST | tr '\n' ' ')

echo ">>> STEP 1: Running samtools mpileup (FULL GENOME)..."
# Filtre yok, pipe yok. Doğrudan diske yazıyoruz.
samtools mpileup -B -Q 20 -f $REF $BAM_FILES > $OUT_DIR/${YEAR}_RAW.mpileup

# Dosya kontrolü
if [ ! -s $OUT_DIR/${YEAR}_RAW.mpileup ]; then
    echo "HATA: mpileup dosyası boş! BAM yollarını veya Referans ID'lerini kontrol et."
    exit 1
fi

echo ">>> STEP 2: Filtering Autosomes (2L, 2R, 3L, 3R, 4)..."
# Ham dosyadan sadece otozomları çekiyoruz. 
# ~ kullanarak olası gizli karakterleri de kapsıyoruz.
awk '$1 ~ /NT_033779\.5|NT_033778\.4|NT_037436\.4|NT_033777\.3|NC_004353\.4/' \
$OUT_DIR/${YEAR}_RAW.mpileup > $OUT_DIR/${YEAR}.mpileup

echo ">>> STEP 3: Renaming chromosomes in the filtered file..."
sed -i 's/NT_033779\.5/2L/g; s/NT_033778\.4/2R/g; s/NT_037436\.4/3L/g; s/NT_033777\.3/3R/g; s/NC_004353\.4/4/g' \
    $OUT_DIR/${YEAR}.mpileup

echo ">>> STEP 4: VarScan mpileup2snp..."
VARSCAN_JAR=$(find $CONDA_PREFIX -name "VarScan*.jar" | head -n 1)
java -Xmx100g -jar $VARSCAN_JAR mpileup2snp $OUT_DIR/${YEAR}.mpileup \
    --min-var-freq 0.01 --p-value 0.05 --min-coverage 20 --min-avg-qual 20 --output-vcf 1 \
    > $OUT_DIR/${YEAR}.vcf

echo ">>> STEP 5: SnpEff Annotation..."
java -Xmx100g -jar $(find $CONDA_PREFIX -name "snpEff.jar" | head -n 1) \
    -v BDGP6.32.105 $OUT_DIR/${YEAR}.vcf > $OUT_DIR/${YEAR}_annotated.vcf

echo ">>> STEP 6: Creating Sync file..."
python $BASE_DIR/vcf_to_sync.py $OUT_DIR/${YEAR}.mpileup $OUT_DIR/${YEAR}_annotated.vcf $OUT_DIR/${YEAR}.sync

# Yer açmak istersen ham dosyayı silebilirsin:
# rm $OUT_DIR/${YEAR}_RAW.mpileup

echo ">>> MEGA ANALİZ TAMAMLANDI: $(date)"


====================================================================
2. STATISTICAL ANALYSIS AND ANCOVA MODELING (R SCRIPT)
====================================================================
# This section details the chronological Fst heatmap generation 
# and the genome-wide ANCOVA model for single nucleotide polymorphisms (SNPs)
# incorporating year, season, and climate variables.


#!/bin/bash
#SBATCH -p orfoz
#SBATCH -A ssenkal
#SBATCH -J mega_popool2_fst
#SBATCH --nodes=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=80G
#SBATCH -t 24:00:00
#SBATCH -o mega_fst_%j.out
#SBATCH -e mega_fst_%j.err

PROGRAMS="/arf/scratch/ssenkal/selin/programs/popoolation2_1201"
OUTDIR="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/fst_results"
SYNC="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/mega_all_years.sync"

mkdir -p $OUTDIR

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate align_env

echo ">>> Popoolation2 Sliding Window Fst Başlıyor..."

perl $PROGRAMS/fst-sliding.pl \
    --input $SYNC \
    --output $OUTDIR/mega_all_years_26_sliding_fst.txt \
    --pool-size 40 \
    --min-count 2 \
    --min-coverage 15 \
    --max-coverage 400 \
    --window-size 5000 \
    --step-size 5000 \
    --suppress-noninformative

echo ">>> İşlem Tamamlandı."


(base) [ssenkal@arf-ui1 fst_results]$ cat PCA_ancova_iklim_fst.R 
library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ============================================================
# 1. CLIMATE VERİSİ
# ============================================================
climate <- data.frame(
  Sample = c("2015_06","2015_08","2015_09","2015_10",
            "2016_06","2016_08","2016_09","2016_10",
            "2017_06","2017_10",
            "2018_06","2018_08","2018_09","2018_10",
            "2019_06","2019_08","2019_09","2019_10",
            "2020_06","2020_08","2020_09","2020_10",
            "2021_06","2021_08","2021_09","2021_10"),
  Month = c("June","August","September","October",
         "June","August","September","October",
         "June","October",
         "June","August","September","October",
         "June","August","September","October",
         "June","August","September","October",
         "June","August","September","October"),
  Year = c(2015,2015,2015,2015,
          2016,2016,2016,2016,
          2017,2017,
          2018,2018,2018,2018,
          2019,2019,2019,2019,
          2020,2020,2020,2020,
          2021,2021,2021,2021),
  MeanTemp  = c(18.2,24.6,23.0,14.1, 21.7,25.7,19.0,13.1, 19.7,11.8,
              21.3,25.4,20.3,14.5, 21.5,24.3,20.2,15.5, 20.2,25.2,23.1,16.6,
              19.1,25.3,18.3,12.5),
  MaxSic  = c(25.4,33.2,31.9,21.2, 30.1,34.3,27.7,21.5, 27.8,19.9,
              30.2,34.2,28.8,22.4, 29.9,33.3,29.1,24.7, 28.5,34.4,32.3,25.3,
              27.0,34.2,26.0,20.7),
  MinSic  = c(12.6,16.7,14.9, 8.6, 13.0,17.2,11.4, 6.4, 12.7, 5.7,
              13.6,17.2,12.3, 8.2, 14.7,15.9,12.1, 8.0, 12.8,16.0,14.8, 9.9,
              12.0,16.3,11.4, 5.7),
  Nem     = c(73.8,47.3,43.4,68.4, 53.4,43.3,50.4,56.7, 63.2,61.7,
              59.6,43.9,48.5,68.7, 69.9,49.9,48.9,60.4, 68.5,42.8,51.1,64.0,
              71.0,47.7,61.3,67.3),
  SicakGun = c(0,28,22,0, 15,27,12,0, 8,0, 17,31,10,0, 15,27,15,2, 12,31,22,3,
               10,30, 3,0),
  SogukGun = c(0,0,0,4, 0,0,4,11, 0,17, 0,0,0,6, 0,0,2,9, 0,0,0,0, 3,0,1,14),
  DonGunu  = c(0,0,0,1, 0,0,0,1, 0,1, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1)
)

# Season etiketi
climate$Season <- ifelse(climate$Month %in% c("June"), "Summer-Start",
                ifelse(climate$Month %in% c("August","September"), "Summer-Summit",
                "Autumn"))

# ============================================================
# 2. FST VERİSİ — sliding window dosyasından genome-wide ortalama
# ============================================================
cat(">>> FST dosyasi okunuyor...\n")
fst_data <- fread("mega_all_years_26_sliding_fst.txt", data.table = FALSE, header = FALSE)

pop_names <- c("2015_06","2015_08","2015_09","2015_10",
               "2016_06","2016_08","2016_09","2016_10",
               "2017_06","2017_10",
               "2018_06","2018_08","2018_09","2018_10",
               "2019_06","2019_08","2019_09","2019_10",
               "2020_06","2020_08","2020_09","2020_10",
               "2021_06","2021_08","2021_09","2021_10")

sample_cols_idx <- 6:ncol(fst_data)
pair_names <- sub("=.*", "", as.character(fst_data[1, sample_cols_idx]))

fst_numeric <- apply(fst_data[, sample_cols_idx], 2, function(col) {
  as.numeric(sub(".*=", "", as.character(col)))
})

mean_fst <- colMeans(fst_numeric, na.rm = TRUE)
names(mean_fst) <- pair_names

# Her örnek için: diğer tüm örneklerle ortalama FST (mean pairwise FST)
fst_long <- data.frame(
  i = as.integer(sub(":.*", "", names(mean_fst))),
  j = as.integer(sub(".*:", "", names(mean_fst))),
  FST = as.numeric(mean_fst)
)

# Her örnek için ortalama pairwise FST
# fst_long'u iki yönlü mirror yaparak her örneği hem i hem j olarak say
fst_mirror2 <- rbind(
  data.frame(Sample = fst_long$i, FST = fst_long$FST),
  data.frame(Sample = fst_long$j, FST = fst_long$FST)
)
mean_pw_fst <- tapply(fst_mirror2$FST, fst_mirror2$Sample, mean, na.rm = TRUE)
climate$MeanFST <- as.numeric(mean_pw_fst[as.character(1:26)])

cat("FST range:", range(climate$MeanFST, na.rm=TRUE), "\n")

# ============================================================
# 3. PCA
# ============================================================
pca_vars <- c("MeanTemp","MaxSic","MinSic","Nem","SicakGun","SogukGun","DonGunu")
pca_input <- climate[, pca_vars]
pca_res <- prcomp(pca_input, scale. = TRUE)

pca_df <- as.data.frame(pca_res$x[, 1:3])
pca_df$Sample  <- climate$Sample
pca_df$Year    <- factor(climate$Year)
pca_df$Season <- climate$Season
pca_df$MeanFST <- climate$MeanFST

var_exp <- round(summary(pca_res)$importance[2, 1:3] * 100, 1)

# PCA plot — yıla göre renk, Seasone göre şekil
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Year, shape = Season, label = Sample)) +
  geom_point(aes(size = MeanFST), alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_size_continuous(range = c(2, 8), name = "Mean pairwise FST") +
  theme_minimal(base_size = 13) +
  labs(
    title = "PCA: CLIMATE Profile (Yesiloz D. melanogaster)",
    subtitle = "Dot Size = mean pairwise FST",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave("PCA_Climate_FST.pdf", p_pca, width = 10, height = 8)
cat(">>> PCA kaydedildi.\n")

# PC yüklemeleri (loadings)
loadings_df <- as.data.frame(pca_res$rotation[, 1:3])
loadings_df$Degisken <- rownames(loadings_df)
cat("\n--- PC Yuklemeleri ---\n")
print(loadings_df)

# ============================================================
# 4. ANCOVA 1: FST ~ Year + CLIMATE değişkenleri
# ============================================================
cat("\n\n============ ANCOVA 1: FST ~ YEAR + CLIMATE ============\n")
ancova_Year <- lm(MeanFST ~ factor(Year) + MeanTemp + Nem + SicakGun + SogukGun,
                 data = climate)
print(summary(ancova_Year))
cat("\n--- ANOVA tablosu ---\n")
print(anova(ancova_Year))

# ============================================================
# 5. ANCOVA 2: FST ~ Season + CLIMATE değişkenleri
# ============================================================
cat("\n\n============ ANCOVA 2: FST ~ SEASON + CLIMATE ============\n")
ancova_Season <- lm(MeanFST ~ Season + MeanTemp + Nem + SicakGun + SogukGun,
                    data = climate)
print(summary(ancova_Season))
cat("\n--- ANOVA tablosu ---\n")
print(anova(ancova_Season))

# ============================================================
# 6. FST ~ PC1 + PC2 regresyon (özet görsel)
# ============================================================
climate$PC1 <- pca_df$PC1
climate$PC2 <- pca_df$PC2

lm_pc <- lm(MeanFST ~ PC1 + PC2, data = climate)
cat("\n\n============ FST ~ PC1 + PC2 ============\n")
print(summary(lm_pc))

p_pc1 <- ggplot(climate, aes(x = PC1, y = MeanFST, color = factor(Year), label = Sample)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40") +
  geom_text_repel(size = 2.8) +
  theme_minimal(base_size = 12) +
  labs(title = "FST ~ PC1 (CLIMATE)", x = paste0("PC1 (", var_exp[1], "%)"),
       y = "Mean pairwise FST", color = "Year")

ggsave("FST_vs_PC1.pdf", p_pc1, width = 9, height = 6)

# ============================================================
# 7. ANCOVA sonuçlarını dosyaya yaz
# ============================================================
sink("ANCOVA_Sonuclari.txt")
cat("========== ANCOVA 1: FST ~ YEAR + CLIMATE ==========\n")
print(summary(ancova_Year))
cat("\nANOVA:\n"); print(anova(ancova_Year))
cat("\n\n========== ANCOVA 2: FST ~ SEASON + CLIMATE ==========\n")
print(summary(ancova_Season))
cat("\nANOVA:\n"); print(anova(ancova_Season))
cat("\n\n========== FST ~ PC1 + PC2 ==========\n")
print(summary(lm_pc))
sink()

cat("\n>>> Tum analizler tamamlandi!\n")
cat("    Ciktilar: PCA_Climate_FST.pdf, FST_vs_PC1.pdf, ANCOVA_Sonuclari.txt\n")



library(data.table)
library(ggplot2)
library(dplyr)

# 1. Pop isimleri (index -> tarih eşlemesi)
pop_names <- c("2015_06", "2015_08", "2015_09", "2015_10",
               "2016_06", "2016_08", "2016_09", "2016_10",
               "2017_06", "2017_10",
               "2018_06", "2018_08", "2018_09", "2018_10",
               "2019_06", "2019_08", "2019_09", "2019_10",
               "2020_06", "2020_08", "2020_09", "2020_10",
               "2021_06", "2021_08", "2021_09", "2021_10")
n_pops <- length(pop_names)

# 2. Load
fst_data <- fread("mega_all_years_26_sliding_fst.txt", data.table = FALSE, header = FALSE)

# 3. Sütun 6'dan itibaren FST
sample_cols_idx <- 6:ncol(fst_data)

# 4. Pair isimlerini ilk satırdan çıkar
pair_names <- sub("=.*", "", as.character(fst_data[1, sample_cols_idx]))

# 5. Sayısal FST değerlerini çıkar
fst_numeric <- apply(fst_data[, sample_cols_idx], 2, function(col) {
  as.numeric(sub(".*=", "", as.character(col)))
})

# 6. Genome-wide ortalama
mean_fst <- colMeans(fst_numeric, na.rm = TRUE)
names(mean_fst) <- pair_names

# 7. Long format - index yerine tarih ismi kullan
fst_long <- data.frame(
  Sample1 = pop_names[as.integer(sub(":.*", "", names(mean_fst)))],
  Sample2 = pop_names[as.integer(sub(".*:", "", names(mean_fst)))],
  FST = as.numeric(mean_fst)
)

# 8. Simetrik mirror
fst_mirror <- data.frame(Sample1 = fst_long$Sample2,
                         Sample2 = fst_long$Sample1,
                         FST = fst_long$FST)
fst_full <- rbind(fst_long, fst_mirror)

# 9. Kronolojik faktör sıralaması
fst_full$Sample1 <- factor(fst_full$Sample1, levels = pop_names)
fst_full$Sample2 <- factor(fst_full$Sample2, levels = pop_names)
fst_full <- fst_full %>% filter(!is.na(FST))

# 10. Plot
ggplot(fst_full, aes(x = Sample1, y = Sample2, fill = FST)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white", high = "red",
    name = expression(italic(F)[ST])
  ) +
  theme_minimal() +
  labs(
    title = "Temporal Pairwise Genetic Differentiation",
    subtitle = "Yesiloz Drosophila melanogaster Population (2015-2021)",
    x = "Sample", y = "Sample"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
    axis.text.y = element_text(size = 8),
    plot.title = element_text(face = "bold", size = 14),
    panel.grid = element_blank()
  )

ggsave("Ankara_Fst_Heatmap_Temporal.pdf", width = 12, height = 10)
cat("Done!\n")


====================================================================
====================================================================

#!/bin/bash
#SBATCH -p hamsi
#SBATCH -A ssenkal
#SBATCH -J SNP_Climate_ANCOVA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH -t 08:00:00
#SBATCH --qos=normal
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/slurm_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/slurm_%j.err

mkdir -p /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

# Install required R packages if missing
Rscript -e "
pkgs <- c('data.table','foreach','doParallel','car')
missing <- pkgs[!pkgs %in% installed.packages()[,'Package']]
if (length(missing) > 0) {
  install.packages(missing, repos='https://cloud.r-project.org')
}
cat('All packages ready.\n')
"

echo ">>> Job started: $(date)"
Rscript /arf/scratch/ssenkal/selin/sekanslar/snp_climate_ancova.R
echo ">>> Job finished: $(date)"


# =============================================================================
# SNP Allele Frequency - Climate ANCOVA Analysis
# Models: Season-only, Year-only, Season*Year, Season+Year+Climate(PCA)
# Input:  mega_all_years_freq_diff_rc
# Output: CLIMATE_ANALYSIS/
# =============================================================================

library(data.table)
library(foreach)
library(doParallel)
library(car)

# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------
registerDoParallel(cores = 40)

INPUT_FILE <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/MAF_ANALYSIS/mega_all_years_freq_diff_rc"
OUT_DIR    <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. SAMPLE METADATA (26 samples, 2017 has only June + October)
# -----------------------------------------------------------------------------
sample_meta <- data.table(
  sample_id = paste0("maa_", 1:26),
  year      = c(rep(2015,4), rep(2016,4), 2017, 2017,
                rep(2018,4), rep(2019,4), rep(2020,4), rep(2021,4)),
  month     = c(6,8,9,10,  6,8,9,10,  6,10,
                6,8,9,10,  6,8,9,10,  6,8,9,10,  6,8,9,10)
)
sample_meta[, season := factor(month, levels=c(6,8,9,10),
                               labels=c("June","August","September","October"))]
sample_meta[, year_val := as.numeric(year)]

cat("Sample metadata:\n")
print(sample_meta)

# -----------------------------------------------------------------------------
# 2. CLIMATE DATA
# Months in dataset: 6=June, 8=August, 9=September, 10=October
# 2017 has no August/September → those rows get NA for climate vars
# -----------------------------------------------------------------------------

# Helper: build a (year x month) lookup table
make_clim <- function(varname, vals_list) {
  rows <- lapply(names(vals_list), function(yr) {
    m <- vals_list[[yr]]
    data.table(year=as.integer(yr), month=as.integer(names(m)),
               value=unlist(m, use.names=FALSE))
  })
  dt <- rbindlist(rows)
  setnames(dt, "value", varname)
  dt
}

avg_temp <- make_clim("avg_temp", list(
  "2015"=c("6"=18.2,"8"=24.6,"9"=23.0,"10"=14.1),
  "2016"=c("6"=21.7,"8"=25.7,"9"=19.0,"10"=13.1),
  "2017"=c("6"=19.7,         "10"=11.8),
  "2018"=c("6"=21.3,"8"=25.4,"9"=20.3,"10"=14.5),
  "2019"=c("6"=21.5,"8"=24.3,"9"=20.2,"10"=15.5),
  "2020"=c("6"=20.2,"8"=25.2,"9"=23.1,"10"=16.6),
  "2021"=c("6"=19.1,"8"=25.3,"9"=18.3,"10"=12.5)))

min_temp <- make_clim("min_temp", list(
  "2015"=c("6"=12.6,"8"=16.7,"9"=14.9,"10"=8.6),
  "2016"=c("6"=13.0,"8"=17.2,"9"=11.4,"10"=6.4),
  "2017"=c("6"=12.7,         "10"=5.7),
  "2018"=c("6"=13.6,"8"=17.2,"9"=12.3,"10"=8.2),
  "2019"=c("6"=14.7,"8"=15.9,"9"=12.1,"10"=8.0),
  "2020"=c("6"=12.8,"8"=16.0,"9"=14.8,"10"=9.9),
  "2021"=c("6"=12.0,"8"=16.3,"9"=11.4,"10"=5.7)))

max_temp <- make_clim("max_temp", list(
  "2015"=c("6"=25.4,"8"=33.2,"9"=31.9,"10"=21.2),
  "2016"=c("6"=30.1,"8"=34.3,"9"=27.7,"10"=21.5),
  "2017"=c("6"=27.8,         "10"=19.9),
  "2018"=c("6"=30.2,"8"=34.2,"9"=28.8,"10"=22.4),
  "2019"=c("6"=29.9,"8"=33.3,"9"=29.1,"10"=24.7),
  "2020"=c("6"=28.5,"8"=34.4,"9"=32.3,"10"=25.3),
  "2021"=c("6"=27.0,"8"=34.2,"9"=26.0,"10"=20.7)))

precip <- make_clim("precipitation", list(
  "2015"=c("6"=152.6,"8"=21.2,"9"=2.7, "10"=21.7),
  "2016"=c("6"=16.1, "8"=4.8, "9"=32.3,"10"=2.4),
  "2017"=c("6"=91.2,          "10"=33.4),
  "2018"=c("6"=22.3, "8"=52.6,"9"=2.0, "10"=67.2),
  "2019"=c("6"=89.6, "8"=9.3, "9"=1.0, "10"=13.8),
  "2020"=c("6"=62.5, "8"=0.4, "9"=8.9, "10"=39.2),
  "2021"=c("6"=118.8,"8"=4.4, "9"=10.4,"10"=25.5)))

hot_days <- make_clim("hot_days", list(
  "2015"=c("6"=0, "8"=28,"9"=22,"10"=0),
  "2016"=c("6"=15,"8"=27,"9"=12,"10"=0),
  "2017"=c("6"=8,         "10"=0),
  "2018"=c("6"=17,"8"=31,"9"=10,"10"=0),
  "2019"=c("6"=15,"8"=27,"9"=15,"10"=2),
  "2020"=c("6"=12,"8"=31,"9"=22,"10"=3),
  "2021"=c("6"=10,"8"=30,"9"=3, "10"=0)))

avg_humidity <- make_clim("avg_humidity", list(
  "2015"=c("6"=73.8,"8"=47.3,"9"=43.4,"10"=68.4),
  "2016"=c("6"=53.4,"8"=43.3,"9"=50.4,"10"=56.7),
  "2017"=c("6"=63.2,         "10"=61.7),
  "2018"=c("6"=59.6,"8"=43.9,"9"=48.5,"10"=68.7),
  "2019"=c("6"=69.9,"8"=49.9,"9"=48.9,"10"=60.4),
  "2020"=c("6"=68.5,"8"=42.8,"9"=51.1,"10"=64.0),
  "2021"=c("6"=71.0,"8"=47.7,"9"=61.3,"10"=67.3)))

min_humidity <- make_clim("min_humidity", list(
  "2015"=c("6"=38.9,"8"=20.4,"9"=18.5,"10"=38.5),
  "2016"=c("6"=22.8,"8"=17.6,"9"=22.6,"10"=29.0),
  "2017"=c("6"=32.5,         "10"=31.7),
  "2018"=c("6"=27.3,"8"=18.5,"9"=23.9,"10"=37.7),
  "2019"=c("6"=34.5,"8"=23.8,"9"=23.3,"10"=31.0),
  "2020"=c("6"=33.7,"8"=20.3,"9"=24.6,"10"=36.3),
  "2021"=c("6"=38.7,"8"=21.8,"9"=35.9,"10"=41.3)))

max_humidity <- make_clim("max_humidity", list(
  "2015"=c("6"=98.9,"8"=77.2,"9"=71.2,"10"=91.9),
  "2016"=c("6"=89.1,"8"=73.2,"9"=80.6,"10"=82.8),
  "2017"=c("6"=89.0,         "10"=85.3),
  "2018"=c("6"=89.8,"8"=72.1,"9"=75.9,"10"=91.4),
  "2019"=c("6"=96.6,"8"=78.0,"9"=76.4,"10"=88.1),
  "2020"=c("6"=94.9,"8"=69.9,"9"=80.6,"10"=86.4),
  "2021"=c("6"=97.7,"8"=76.0,"9"=87.2,"10"=91.4)))

# Merge all climate tables
clim_all <- Reduce(function(a,b) merge(a,b,by=c("year","month"),all=TRUE),
                   list(avg_temp, min_temp, max_temp, precip,
                        hot_days, avg_humidity, min_humidity, max_humidity))

# Derived variables
clim_all[, temp_amplitude   := max_temp - min_temp]
clim_all[, humidity_range   := max_humidity - min_humidity]

# Year-to-year delta (previous year same month)
setorder(clim_all, month, year)
clim_all[, delta_avg_temp   := avg_temp   - shift(avg_temp),   by=month]
clim_all[, delta_min_temp   := min_temp   - shift(min_temp),   by=month]
clim_all[, delta_max_temp   := max_temp   - shift(max_temp),   by=month]
clim_all[, delta_humidity   := avg_humidity - shift(avg_humidity), by=month]
clim_all[, delta_precip     := precipitation - shift(precipitation), by=month]

CLIM_VARS <- c("avg_temp","min_temp","max_temp","precipitation",
               "hot_days","avg_humidity","min_humidity","max_humidity",
               "temp_amplitude","humidity_range",
               "delta_avg_temp","delta_min_temp","delta_max_temp",
               "delta_humidity","delta_precip")

# Merge climate into sample metadata
sample_meta <- merge(sample_meta, clim_all[, c("year","month",CLIM_VARS), with=FALSE],
                     by=c("year","month"), all.x=TRUE)

# -----------------------------------------------------------------------------
# 3. PCA ON CLIMATE VARIABLES (across the 26 samples)
# -----------------------------------------------------------------------------
clim_matrix <- as.matrix(sample_meta[, CLIM_VARS, with=FALSE])
# Use only complete rows for PCA, then predict for all
complete_rows <- complete.cases(clim_matrix)
pca_fit <- prcomp(clim_matrix[complete_rows,], scale.=TRUE, center=TRUE)

cat("\nPCA variance explained:\n")
pct <- summary(pca_fit)$importance[2,] * 100
print(round(pct[1:5], 1))

# Predict PC scores for all 26 samples (NA rows get NA scores)
pc_scores <- predict(pca_fit, newdata=clim_matrix)
sample_meta[, PC1 := pc_scores[,1]]
sample_meta[, PC2 := pc_scores[,2]]
sample_meta[, PC3 := pc_scores[,3]]

cat("\nPC1+PC2 variance explained:", round(sum(pct[1:2]),1), "%\n")

# -----------------------------------------------------------------------------
# 4. READ INPUT FILE
# -----------------------------------------------------------------------------
cat("\n>>> Reading input file...\n")

# Header line starts with ##; fread skips it as a comment.
# Read column names manually, then pass directly to fread.
raw_header <- readLines(INPUT_FILE, n = 1)
col_names  <- strsplit(sub("^##", "", raw_header), "\t")[[1]]
cat("Columns detected:", length(col_names), "\n")

# skip=1 skips the ## header line; we supply col.names ourselves
dt <- fread(INPUT_FILE, skip = 1, header = FALSE, col.names = col_names)

# Parse allele frequencies for each sample
cat(">>> Parsing allele frequencies...\n")
y_list <- vector("list", nrow(sample_meta))

for (i in seq_len(nrow(sample_meta))) {
  samp  <- sample_meta$sample_id[i]
  if (!samp %in% names(dt)) next
  counts <- tstrsplit(dt[[samp]], "/", type.convert=TRUE)
  tot    <- counts[[2]]
  af     <- round(counts[[1]] / tot, 6)
  
  tmp <- dt[, .(chr, pos)]
  tmp[, AF       := af]
  tmp[, tot_c    := tot]
  tmp[, season   := sample_meta$season[i]]
  tmp[, year_val := sample_meta$year_val[i]]
  tmp[, PC1      := sample_meta$PC1[i]]
  tmp[, PC2      := sample_meta$PC2[i]]
  tmp[, PC3      := sample_meta$PC3[i]]
  tmp[, avg_temp     := sample_meta$avg_temp[i]]
  tmp[, avg_humidity := sample_meta$avg_humidity[i]]
  tmp[, temp_amplitude := sample_meta$temp_amplitude[i]]
  
  y_list[[i]] <- tmp[tot_c >= 15]
  rm(tmp, counts, af, tot); gc()
}

rm(dt); gc()
y_master <- rbindlist(y_list)
rm(y_list); gc()

snps <- unique(y_master[, .(chr, pos)])
cat("Total SNPs:", nrow(snps), "\n")

# -----------------------------------------------------------------------------
# 5. PARALLEL MODEL FITTING
# -----------------------------------------------------------------------------
cat("\n>>> Running parallel ANCOVA models...\n")

res <- foreach(i = 1:nrow(snps),
               .combine  = rbind,
               .packages = c("data.table","car"),
               .errorhandling = "remove") %dopar% {

  s   <- snps[i]
  sub <- y_master[chr == s$chr & pos == s$pos]
  if (nrow(sub) < 10) return(NULL)

  tryCatch({
    # ── Model 1: Season only ──────────────────────────────────────────────
    m1   <- lm(AF ~ season, data=sub)
    p1   <- anova(m1)$`Pr(>F)`[1]

    # ── Model 2: Year only ───────────────────────────────────────────────
    m2   <- lm(AF ~ year_val, data=sub)
    p2   <- anova(m2)$`Pr(>F)`[1]

    # ── Model 3: Season * Year (ANCOVA) ──────────────────────────────────
    m3   <- lm(AF ~ season * year_val, data=sub)
    a3   <- Anova(m3, type="III")
    p3_season <- a3$`Pr(>F)`[2]   # season
    p3_year   <- a3$`Pr(>F)`[3]   # year_val
    p3_inter  <- a3$`Pr(>F)`[4]   # season:year_val

    # ── Model 4: Season + Year + Climate PCA ─────────────────────────────
    # Use PC1+PC2; drop rows where PCs are NA (2017 Aug/Sep)
    sub4 <- sub[!is.na(PC1) & !is.na(PC2)]
    if (nrow(sub4) < 10) {
      p4_season <- NA_real_; p4_year <- NA_real_
      p4_PC1 <- NA_real_;    p4_PC2 <- NA_real_
    } else {
      m4   <- lm(AF ~ season + year_val + PC1 + PC2, data=sub4)
      a4   <- Anova(m4, type="III")
      p4_season <- a4$`Pr(>F)`[2]
      p4_year   <- a4$`Pr(>F)`[3]
      p4_PC1    <- a4$`Pr(>F)`[4]
      p4_PC2    <- a4$`Pr(>F)`[5]
    }

    data.table(
      chr=s$chr, pos=s$pos,
      p1_season=p1,
      p2_year=p2,
      p3_season=p3_season, p3_year=p3_year, p3_interaction=p3_inter,
      p4_season=p4_season, p4_year=p4_year,
      p4_PC1=p4_PC1, p4_PC2=p4_PC2
    )
  }, error=function(e) NULL)
}

cat("Models done. Applying BH correction...\n")

# -----------------------------------------------------------------------------
# 6. BH CORRECTION (separately per test column)
# -----------------------------------------------------------------------------
p_cols <- c("p1_season","p2_year",
            "p3_season","p3_year","p3_interaction",
            "p4_season","p4_year","p4_PC1","p4_PC2")

for (pc in p_cols) {
  qc <- sub("^p", "q", pc)
  res[, (qc) := p.adjust(get(pc), method="BH")]
}

# -----------------------------------------------------------------------------
# 7. SAVE FULL RESULTS
# -----------------------------------------------------------------------------
full_out <- file.path(OUT_DIR, "snp_climate_ancova_full.tsv")
fwrite(res, full_out, sep="\t")
cat("Full results saved:", full_out, "\n")

# -----------------------------------------------------------------------------
# 8. SIGNIFICANT SNP SUBSETS (q < 0.05 per test)
# -----------------------------------------------------------------------------
sig_tests <- list(
  season_only      = res[q1_season      < 0.05, .(chr, pos, p1_season, q1_season)],
  year_only        = res[q2_year        < 0.05, .(chr, pos, p2_year,   q2_year)],
  ancova_season    = res[q3_season      < 0.05, .(chr, pos, p3_season, q3_season,
                                                   p3_year, q3_year,
                                                   p3_interaction, q3_interaction)],
  ancova_year      = res[q3_year        < 0.05, .(chr, pos, p3_season, q3_season,
                                                   p3_year, q3_year,
                                                   p3_interaction, q3_interaction)],
  ancova_interact  = res[q3_interaction < 0.05, .(chr, pos, p3_season, q3_season,
                                                   p3_year, q3_year,
                                                   p3_interaction, q3_interaction)],
  climate_PC1      = res[q4_PC1         < 0.05, .(chr, pos, p4_season, q4_season,
                                                   p4_year, q4_year,
                                                   p4_PC1, q4_PC1,
                                                   p4_PC2, q4_PC2)],
  climate_PC2      = res[q4_PC2         < 0.05, .(chr, pos, p4_season, q4_season,
                                                   p4_year, q4_year,
                                                   p4_PC1, q4_PC1,
                                                   p4_PC2, q4_PC2)]
)

for (nm in names(sig_tests)) {
  out_path <- file.path(OUT_DIR, paste0("sig_", nm, ".tsv"))
  fwrite(sig_tests[[nm]], out_path, sep="\t")
  cat(sprintf("  %-20s → %d significant SNPs → %s\n",
              nm, nrow(sig_tests[[nm]]), out_path))
}

# -----------------------------------------------------------------------------
# 9. PCA LOADINGS REPORT
# -----------------------------------------------------------------------------
loadings_out <- file.path(OUT_DIR, "climate_pca_loadings.tsv")
loadings_dt  <- as.data.table(pca_fit$rotation[, 1:3], keep.rownames=TRUE)
setnames(loadings_dt, "rn", "climate_variable")
fwrite(loadings_dt, loadings_out, sep="\t")
cat("PCA loadings saved:", loadings_out, "\n")

cat("\n>>> All done:", format(Sys.time()), "\n")



====================================================================
3. FUNCTIONAL GO ENRICHMENT ANALYSIS (PYTHON SCRIPT)
====================================================================
# Python implementation of gseapy for Gene Ontology (GO) enrichment
# utilizing the year-significant SNP gene list dataset.


GO Enrichment Analysis using gseapy (offline/local gene sets)
Input:  year_sig_gene_list.txt
Output: GO_ENRICHMENT/ directory
"""

import gseapy as gp
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

GENE_FILE = "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/year_sig_gene_list.txt"
OUT_DIR   = "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/GO_ENRICHMENT"
os.makedirs(OUT_DIR, exist_ok=True)

# Read gene list
with open(GENE_FILE) as f:
    genes = [line.strip() for line in f if line.strip()]
print(f"Input genes: {len(genes)}")

# Local gene sets for Fly
gene_sets = {
    "BP": "GO_Biological_Process_2018",
    "MF": "GO_Molecular_Function_2018",
    "CC": "GO_Cellular_Component_2018",
    "KEGG": "KEGG_2019",
}

all_results = []

for cat, gs in gene_sets.items():
    print(f"\n>>> Running enrichment: {gs}")
    try:
        enr = gp.enrichr(
            gene_list     = genes,
            gene_sets     = gs,
            organism      = "fly",
            outdir        = None,
            cutoff        = 0.05,
            no_plot       = True,
        )

        if enr.results is not None and len(enr.results) > 0:
            sig = enr.results[enr.results["Adjusted P-value"] < 0.05].copy()
            sig["Ontology"] = cat
            print(f"  Significant terms: {len(sig)}")

            if len(sig) > 0:
                all_results.append(sig)
                sig.to_csv(f"{OUT_DIR}/GO_{cat}_results.tsv", sep="\t", index=False)

                # Dot plot
                top20 = sig.head(20).copy()
                top20["gene_count"] = top20["Overlap"].apply(lambda x: int(x.split("/")[0]))
                top20["-log10(padj)"] = -np.log10(top20["Adjusted P-value"].clip(lower=1e-300))

                fig, ax = plt.subplots(figsize=(11, max(6, len(top20)*0.45)))
                scatter = ax.scatter(
                    top20["-log10(padj)"],
                    range(len(top20)),
                    s=top20["gene_count"] * 25,
                    c=top20["Adjusted P-value"],
                    cmap="RdYlGn_r",
                    alpha=0.85,
                    edgecolors="grey",
                    linewidth=0.5
                )
                ax.set_yticks(range(len(top20)))
                ax.set_yticklabels(top20["Term"], fontsize=9)
                ax.set_xlabel("-log10(Adjusted P-value)", fontsize=11)
                ax.set_title(f"GO {cat} Enrichment\n(Top {len(top20)} significant terms, q<0.05)",
                             fontsize=12, fontweight="bold")
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label("Adjusted P-value", fontsize=9)
                plt.tight_layout()
                plt.savefig(f"{OUT_DIR}/dotplot_{cat}.pdf", bbox_inches="tight")
                plt.savefig(f"{OUT_DIR}/dotplot_{cat}.png", dpi=150, bbox_inches="tight")
                plt.close()

                # Bar plot
                fig, ax = plt.subplots(figsize=(11, max(6, len(top20)*0.45)))
                colors = plt.cm.RdYlGn_r(
                    (top20["Adjusted P-value"] - top20["Adjusted P-value"].min()) /
                    (top20["Adjusted P-value"].max() - top20["Adjusted P-value"].min() + 1e-10)
                )
                bars = ax.barh(range(len(top20)), top20["-log10(padj)"],
                               color=colors, edgecolor="grey", linewidth=0.5)
                ax.set_yticks(range(len(top20)))
                ax.set_yticklabels(top20["Term"], fontsize=9)
                ax.set_xlabel("-log10(Adjusted P-value)", fontsize=11)
                ax.set_title(f"GO {cat} Enrichment\n(Top {len(top20)} significant terms, q<0.05)",
                             fontsize=12, fontweight="bold")
                ax.invert_yaxis()
                plt.tight_layout()
                plt.savefig(f"{OUT_DIR}/barplot_{cat}.pdf", bbox_inches="tight")
                plt.savefig(f"{OUT_DIR}/barplot_{cat}.png", dpi=150, bbox_inches="tight")
                plt.close()

                print(f"  Plots saved: dotplot_{cat}.pdf, barplot_{cat}.pdf")
        else:
            print(f"  No significant terms")

    except Exception as e:
        print(f"  Error: {e}")

# Combined
if all_results:
    combined = pd.concat(all_results, ignore_index=True)
    combined = combined.sort_values("Adjusted P-value")
    combined.to_csv(f"{OUT_DIR}/GO_all_significant.tsv", sep="\t", index=False)
    print(f"\nTotal significant GO terms: {len(combined)}")

    # Summary bubble plot - all categories
    top10_each = combined.groupby("Ontology").head(10).copy()
    top10_each["gene_count"] = top10_each["Overlap"].apply(lambda x: int(x.split("/")[0]))
    top10_each["-log10(padj)"] = -np.log10(top10_each["Adjusted P-value"].clip(lower=1e-300))

    color_map = {"BP": "#E74C3C", "MF": "#3498DB", "CC": "#2ECC71", "KEGG": "#9B59B6"}
    fig, ax = plt.subplots(figsize=(13, max(8, len(top10_each)*0.4)))
    for _, row in top10_each.iterrows():
        ax.scatter(row["-log10(padj)"], row["Term"],
                   s=row["gene_count"] * 30,
                   c=color_map.get(row["Ontology"], "grey"),
                   alpha=0.8, edgecolors="grey", linewidth=0.5)
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0],[0], marker='o', color='w',
                               markerfacecolor=v, markersize=10, label=k)
                       for k, v in color_map.items() if k in top10_each["Ontology"].values]
    ax.legend(handles=legend_elements, title="Ontology", fontsize=9)
    ax.set_xlabel("-log10(Adjusted P-value)", fontsize=11)
    ax.set_title("GO Enrichment Summary — Year-Significant SNP Genes\n(Top 10 per category, q<0.05)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/GO_summary_bubble.pdf", bbox_inches="tight")
    plt.savefig(f"{OUT_DIR}/GO_summary_bubble.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Summary bubble plot saved.")
else:
    print("\nNo significant terms found.")

print(f"\n>>> Done: all results in {OUT_DIR}")
