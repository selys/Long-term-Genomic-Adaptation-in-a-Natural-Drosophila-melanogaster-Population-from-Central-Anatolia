#!/bin/bash
#SBATCH -p hamsi
#SBATCH -A ssenkal
#SBATCH -J gene_length_HM
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=32G
#SBATCH -t 02:00:00
#SBATCH --qos=normal
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/gene_length_HM_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/gene_length_HM_%j.err

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

BASE_OUT="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL"
OUT_DIR="$BASE_OUT/HIGH_MODERATE_ONLY"
mkdir -p $OUT_DIR

# -----------------------------------------------------------------------
# Reuse everything already produced by gene_length_control.sh -- no need
# to re-download the GFF or re-annotate the VCF.
# -----------------------------------------------------------------------
GENES_GTF="$BASE_OUT/genes_raw.gtf"
FULL_VCF="$BASE_OUT/mega_all_years_reannotated.vcf"
SIG_VCF="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/year_sig_annotated.vcf.gz"

if [ ! -f "$GENES_GTF" ]; then
    echo "HATA: $GENES_GTF yok. Once gene_length_control.sh'i calistir."
    exit 1
fi
if [ ! -f "$FULL_VCF" ]; then
    echo "HATA: $FULL_VCF yok. Once gene_length_control.sh'i calistir (reannotated VCF'i uretir)."
    exit 1
fi
echo "Using GENES_GTF: $GENES_GTF"
echo "Using FULL_VCF:  $FULL_VCF"
echo "Using SIG_VCF:   $SIG_VCF"

# -----------------------------------------------------------------------
# Extract SNP-gene pairs, but ONLY where Annotation_Impact (ANN field 3)
# is HIGH or MODERATE -- this is the exact SNP universe that feeds GO
# enrichment (the 500-variant / 395-gene HIGH/MODERATE set), not all
# 9,055 year-significant SNPs regardless of impact.
# -----------------------------------------------------------------------
AUTOS="2L 2R 3L 3R 4"

echo ">>> Extracting HIGH/MODERATE background SNP-gene pairs..."
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' "$FULL_VCF" | \
  awk -F'\t' -v autos="$AUTOS" '
    BEGIN { split(autos, a, " "); for (i in a) keep[a[i]] = 1 }
    {
      if (!($1 in keep)) next;
      n=split($3, anns, ",");
      for(i=1;i<=n;i++){
        split(anns[i], f, "|");
        impact=f[3]; gene=f[4];
        if((impact=="HIGH" || impact=="MODERATE") && gene!="") print $1"\t"$2"\t"gene;
      }
    }' | sort -u > $OUT_DIR/background_HM_chr_pos_gene.txt

echo ">>> Extracting HIGH/MODERATE year-significant SNP-gene pairs..."
bcftools query -f '%CHROM\t%POS\t%INFO/ANN\n' "$SIG_VCF" | \
  awk -F'\t' -v autos="$AUTOS" '
    BEGIN { split(autos, a, " "); for (i in a) keep[a[i]] = 1 }
    {
      if (!($1 in keep)) next;
      n=split($3, anns, ",");
      for(i=1;i<=n;i++){
        split(anns[i], f, "|");
        impact=f[3]; gene=f[4];
        if((impact=="HIGH" || impact=="MODERATE") && gene!="") print $1"\t"$2"\t"gene;
      }
    }' | sort -u > $OUT_DIR/sig_HM_chr_pos_gene.txt

echo "Background HIGH/MODERATE unique chr+pos+gene rows: $(wc -l < $OUT_DIR/background_HM_chr_pos_gene.txt)"
echo "Year-sig HIGH/MODERATE unique chr+pos+gene rows:   $(wc -l < $OUT_DIR/sig_HM_chr_pos_gene.txt)"

if [ ! -s "$OUT_DIR/background_HM_chr_pos_gene.txt" ]; then
    echo "HATA: background_HM_chr_pos_gene.txt bos -- ANN alaninda HIGH/MODERATE bulunamadi mi kontrol et."
    exit 1
fi

# -----------------------------------------------------------------------
# R: repeat the length-bias tests restricted to this HIGH/MODERATE universe,
# plus a direct test of whether the 395 candidate genes are themselves
# systematically longer than the genome background.
# -----------------------------------------------------------------------
cat > $OUT_DIR/gene_length_HM_analysis.R << 'RSCRIPT'
library(data.table)
library(ggplot2)

BASE_OUT <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL"
OUT_DIR  <- file.path(BASE_OUT, "HIGH_MODERATE_ONLY")

# --- Gene lengths (reuse the already-parsed, autosome-filtered GTF logic) ---
gtf <- fread(file.path(BASE_OUT, "genes_raw.gtf"), header = FALSE, sep = "\t", quote = "")
setnames(gtf, c("chr","source","feature","start","end","score","strand","frame","attr"))
gtf[, length_bp := end - start + 1]

extract_attr <- function(attr, key) {
  gtf_pat <- paste0(key, ' "([^"]+)"')
  gff_pat <- paste0('(^|;)\\s*', key, '=([^;]+)')
  out <- rep(NA_character_, length(attr))
  hit_gtf <- grepl(gtf_pat, attr)
  out[hit_gtf] <- sub(paste0('.*', gtf_pat, '.*'), '\\1', attr[hit_gtf])
  hit_gff <- is.na(out) & grepl(gff_pat, attr)
  out[hit_gff] <- sub(paste0('.*', gff_pat, '.*'), '\\2', attr[hit_gff])
  out
}
gtf[, gene_name := extract_attr(attr, "gene_name")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "gene")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "Name")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "gene_id")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "ID")]
gtf <- gtf[!is.na(gene_name)]

AUTOSOME_ACCESSIONS <- c("NT_033779.5","NT_033778.4","NT_037436.4","NT_033777.3","NC_004353.4")
gtf <- gtf[chr %in% AUTOSOME_ACCESSIONS]
gene_lengths <- gtf[, .(length_bp = max(length_bp)), by = gene_name]
cat("Autosomal genes with length info:", nrow(gene_lengths), "\n")

# --- HIGH/MODERATE background & significant SNP-gene tables ---
bg  <- fread(file.path(OUT_DIR, "background_HM_chr_pos_gene.txt"), header = FALSE,
             col.names = c("chr","pos","gene_name"))
sig <- fread(file.path(OUT_DIR, "sig_HM_chr_pos_gene.txt"), header = FALSE,
             col.names = c("chr","pos","gene_name"))

bg_counts  <- bg[,  .(n_tested_HM = .N), by = gene_name]
sig_counts <- sig[, .(n_sig_HM    = .N), by = gene_name]

merged <- merge(bg_counts, sig_counts, by = "gene_name", all.x = TRUE)
merged[is.na(n_sig_HM), n_sig_HM := 0]
merged <- merge(merged, gene_lengths, by = "gene_name")
merged <- merged[length_bp > 0 & n_tested_HM > 0]
merged[, prop_sig_HM := n_sig_HM / n_tested_HM]
merged[, density_bp_HM := n_sig_HM / length_bp]
cat("Genes in HIGH/MODERATE merged table:", nrow(merged), "\n")

# --- Same 3 tests, restricted to the HIGH/MODERATE universe ---
ct1 <- cor.test(merged$length_bp, merged$n_sig_HM, method = "spearman")
ct2 <- cor.test(merged$length_bp, merged$prop_sig_HM, method = "spearman")
merged[, log_length := log(length_bp)]
pois_fit <- glm(n_sig_HM ~ log_length + offset(log(n_tested_HM)), data = merged, family = poisson())

# --- Direct test: are the 395 candidate (year-sig HIGH/MODERATE) genes
#     themselves systematically LONGER than the rest of the autosomal
#     gene universe? (Wilcoxon rank-sum / Mann-Whitney U test)
candidate_genes <- unique(sig$gene_name)
gene_lengths[, is_candidate := gene_name %in% candidate_genes]
cat("\nCandidate (HIGH/MODERATE year-sig) genes found in length table:",
    sum(gene_lengths$is_candidate), "of", length(candidate_genes), "\n")

wt <- wilcox.test(length_bp ~ is_candidate, data = gene_lengths)
med_candidate  <- median(gene_lengths[is_candidate == TRUE, length_bp])
med_background <- median(gene_lengths[is_candidate == FALSE, length_bp])

sink(file.path(OUT_DIR, "gene_length_HM_summary.txt"))
cat("Gene-length confound check -- HIGH/MODERATE (GO-enrichment) SNP set only\n")
cat("==========================================================================\n\n")
cat("Genes analyzed (HIGH/MODERATE universe):", nrow(merged), "\n\n")

cat("Test 1 - raw HIGH/MODERATE count vs length (Spearman rho, p):\n")
cat(sprintf("  rho = %.4f, p = %.3g\n\n", ct1$estimate, ct1$p.value))

cat("Test 2 - proportion significant (HIGH/MODERATE) vs length (Spearman rho, p):\n")
cat(sprintf("  rho = %.4f, p = %.3g\n\n", ct2$estimate, ct2$p.value))

cat("Test 3 - Poisson regression coefficient on log(length), offset=log(n_tested_HM):\n")
print(coef(summary(pois_fit)))

cat("\n--------------------------------------------------------------------\n")
cat("Direct test: are the", length(candidate_genes), "candidate genes (year-sig\n")
cat("HIGH/MODERATE) systematically LONGER than the rest of the autosomal genome?\n")
cat(sprintf("  Median length, candidate genes:   %.0f bp\n", med_candidate))
cat(sprintf("  Median length, background genes:  %.0f bp\n", med_background))
cat("  Wilcoxon rank-sum test:\n")
print(wt)
sink()

fwrite(merged, file.path(OUT_DIR, "gene_length_vs_significance_HM.tsv"), sep = "\t")
fwrite(gene_lengths[, .(gene_name, length_bp, is_candidate)],
       file.path(OUT_DIR, "candidate_vs_background_length.tsv"), sep = "\t")

p1 <- ggplot(gene_lengths, aes(x = length_bp, fill = is_candidate)) +
  geom_density(alpha = 0.5) +
  scale_x_log10() +
  labs(title = "Gene length distribution: HIGH/MODERATE candidate genes vs background",
       x = "Gene length (bp, log10 scale)", y = "Density",
       fill = "Is candidate gene") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT_DIR, "candidate_vs_background_length_density.pdf"), p1, width = 7, height = 5)

cat("\n>>> Done. Outputs in", OUT_DIR, "\n")
RSCRIPT

Rscript $OUT_DIR/gene_length_HM_analysis.R

echo ">>> Job finished: $(date)"
