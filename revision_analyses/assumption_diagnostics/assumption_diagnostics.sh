#!/bin/bash
#SBATCH -p orfoz
#SBATCH -A ssenkal
#SBATCH -J assumption_diagnostics
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=48G
#SBATCH -t 06:00:00
#SBATCH --qos=normal
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/assumptions_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/assumptions_%j.err

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

mkdir -p /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/ASSUMPTION_CHECK

cat > /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/ASSUMPTION_CHECK/assumption_check.R << 'RSCRIPT'
# =============================================================================
# Parametric assumption diagnostics + non-parametric robustness check
# Addresses Reviewer 2: "not clear if the assumptions of the parametric tests
# are fulfilled (as both AF and FST range between 0 and 1)"
# =============================================================================
library(data.table)
library(foreach)
library(doParallel)
library(ggplot2)

set.seed(42)
registerDoParallel(cores = 16)

INPUT_FILE <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/MAF_ANALYSIS/mega_all_years_freq_diff_rc"
SIG_FILE   <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/sig_year_only.tsv"
OUT_DIR    <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/ASSUMPTION_CHECK"

# -----------------------------------------------------------------------
# 1. Sample metadata (same as main ANCOVA script)
# -----------------------------------------------------------------------
sample_meta <- data.table(
  sample_id = paste0("maa_", 1:26),
  year  = c(rep(2015,4), rep(2016,4), 2017, 2017, rep(2018,4), rep(2019,4), rep(2020,4), rep(2021,4)),
  month = c(6,8,9,10, 6,8,9,10, 6,10, 6,8,9,10, 6,8,9,10, 6,8,9,10, 6,8,9,10)
)
sample_meta[, season := factor(month, levels=c(6,8,9,10), labels=c("June","August","September","October"))]
sample_meta[, year_val := as.numeric(year)]

# -----------------------------------------------------------------------
# 2. Read input, take a random sample of SNPs (2000) + all year-significant
#    SNPs (or a 2000-SNP subsample of them if the list is very large), so
#    diagnostics cover both the "typical" SNP and the "significant" SNP case
# -----------------------------------------------------------------------
raw_header <- readLines(INPUT_FILE, n = 1)
col_names <- strsplit(sub("^##", "", raw_header), "\t")[[1]]
dt <- fread(INPUT_FILE, skip = 1, header = FALSE, col.names = col_names)
cat("Total SNPs in file:", nrow(dt), "\n")

set.seed(42)
random_idx <- sample(seq_len(nrow(dt)), min(2000, nrow(dt)))

sig <- fread(SIG_FILE)
sig_key <- paste(sig$chr, sig$pos)
dt_key  <- paste(dt$chr, dt$pos)
sig_idx <- which(dt_key %in% sig_key)
if (length(sig_idx) > 2000) sig_idx <- sample(sig_idx, 2000)

# Build each group as its own tagged table FIRST, then resolve overlap
# explicitly. Do not rely on row position after combining/unique() --
# a SNP that lands in both the random draw and the significant set must
# not silently end up mislabeled as "random". We tag it "year_sig" when
# it belongs to both, since that is the more informative label.
random_rows <- dt[random_idx]
random_rows[, snp_group := "random"]

sig_rows <- dt[sig_idx]
sig_rows[, snp_group := "year_sig"]

combined <- rbindlist(list(sig_rows, random_rows))  # sig_rows first = wins ties
dt_sub <- unique(combined, by = c("chr", "pos"))     # one row per chr+pos, correct label

n_overlap <- length(random_idx) + length(sig_idx) - nrow(dt_sub)
cat("Diagnostic subsample size:", nrow(dt_sub),
    " (random:", length(random_idx), ", year_sig:", length(sig_idx),
    ", overlap re-labeled as year_sig:", n_overlap, ")\n")

# -----------------------------------------------------------------------
# 3. Build long-format AF table for the subsample
# -----------------------------------------------------------------------
y_list <- vector("list", nrow(sample_meta))
for (i in seq_len(nrow(sample_meta))) {
  samp <- sample_meta$sample_id[i]
  counts <- tstrsplit(dt_sub[[samp]], "/", type.convert = TRUE)
  tot <- counts[[2]]
  af  <- round(counts[[1]] / tot, 6)
  tmp <- dt_sub[, .(chr, pos, snp_group)]
  tmp[, AF := af]
  tmp[, tot_c := tot]
  tmp[, season := sample_meta$season[i]]
  tmp[, year_val := sample_meta$year_val[i]]
  y_list[[i]] <- tmp[tot_c >= 15]
}
y_master <- rbindlist(y_list)
snps <- unique(y_master[, .(chr, pos, snp_group)])
cat("SNPs with usable coverage in subsample:", nrow(snps), "\n")

# -----------------------------------------------------------------------
# 4. Per-SNP diagnostics: refit Model 2 (AF ~ Year), collect
#    (a) Shapiro-Wilk p-value on residuals
#    (b) parametric p-value (lm) vs non-parametric p-value (Kruskal-Wallis
#        on AF by year, and Spearman AF~year_val) for concordance check
# -----------------------------------------------------------------------
diag_res <- foreach(i = 1:nrow(snps), .combine = rbind,
                     .packages = "data.table", .errorhandling = "remove") %dopar% {
  s <- snps[i]
  sub <- y_master[chr == s$chr & pos == s$pos]
  if (nrow(sub) < 10) return(NULL)

  m2 <- lm(AF ~ year_val, data = sub)
  resid_m2 <- residuals(m2)
  sw <- tryCatch(shapiro.test(resid_m2)$p.value, error = function(e) NA_real_)
  p_lm <- anova(m2)$`Pr(>F)`[1]

  kw  <- tryCatch(kruskal.test(AF ~ factor(year_val), data = sub)$p.value, error = function(e) NA_real_)
  sp  <- tryCatch(cor.test(sub$AF, sub$year_val, method = "spearman")$p.value, error = function(e) NA_real_)

  data.table(chr = s$chr, pos = s$pos, snp_group = s$snp_group,
             shapiro_p = sw, p_lm_year = p_lm, p_kruskal_year = kw, p_spearman_year = sp)
}

diag_res[, q_lm_year := p.adjust(p_lm_year, method = "BH")]
diag_res[, q_kruskal_year := p.adjust(p_kruskal_year, method = "BH")]
diag_res[, q_spearman_year := p.adjust(p_spearman_year, method = "BH")]

fwrite(diag_res, file.path(OUT_DIR, "assumption_diagnostics_full.tsv"), sep = "\t")

# -----------------------------------------------------------------------
# 5. Summaries
# -----------------------------------------------------------------------
sink(file.path(OUT_DIR, "assumption_check_summary.txt"))
cat("Parametric assumption diagnostics\n==================================\n\n")

cat("Residual normality (Shapiro-Wilk on lm(AF ~ Year) residuals):\n")
cat(sprintf("  %% of SNPs with Shapiro p < 0.05 (non-normal residuals): %.1f%%\n\n",
            100 * mean(diag_res$shapiro_p < 0.05, na.rm = TRUE)))

cat("Concordance between parametric (lm) and non-parametric tests, at q<0.05:\n")
lm_sig  <- diag_res$q_lm_year < 0.05
kw_sig  <- diag_res$q_kruskal_year < 0.05
sp_sig  <- diag_res$q_spearman_year < 0.05

cat(sprintf("  lm-significant SNPs:        %d / %d\n", sum(lm_sig, na.rm=TRUE), nrow(diag_res)))
cat(sprintf("  Kruskal-significant SNPs:   %d / %d\n", sum(kw_sig, na.rm=TRUE), nrow(diag_res)))
cat(sprintf("  Spearman-significant SNPs:  %d / %d\n", sum(sp_sig, na.rm=TRUE), nrow(diag_res)))
cat(sprintf("  Agreement lm vs Kruskal (both sig or both non-sig): %.1f%%\n",
            100 * mean(lm_sig == kw_sig, na.rm = TRUE)))
cat(sprintf("  Agreement lm vs Spearman (both sig or both non-sig): %.1f%%\n",
            100 * mean(lm_sig == sp_sig, na.rm = TRUE)))

cat("\nBreakdown by SNP group (random background vs year-significant):\n")
print(diag_res[, .(n = .N,
                    pct_nonnormal_resid = round(100*mean(shapiro_p < 0.05, na.rm=TRUE),1),
                    pct_lm_sig = round(100*mean(q_lm_year < 0.05, na.rm=TRUE),1),
                    pct_kruskal_sig = round(100*mean(q_kruskal_year < 0.05, na.rm=TRUE),1)),
                by = snp_group])
sink()

# -----------------------------------------------------------------------
# 6. FST ~ PC1 + PC2 residual diagnostics (n=22 cohort-level model)
#    Rebuilds the SAME 15-variable delta-based climate PCA used in
#    snp_climate_ancova.R (the script that actually produced Fig. 2's
#    49.2%/19.6% PCA and n=22 regression), so PC1/PC2 here are identical
#    to the manuscript's, not re-derived from a different/older script.
# -----------------------------------------------------------------------
FST_FILE <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/fst_results/mega_all_years_25_sliding_fst_no2015_09.txt"

if (!file.exists(FST_FILE)) {
  cat("\n>>> FST_FILE not found, skipping FST~PC1+PC2 diagnostics:", FST_FILE, "\n")
} else {

  cat("\n>>> Rebuilding 15-variable climate PCA (same as snp_climate_ancova.R) ...\n")

  make_clim <- function(varname, vals_list) {
    rows <- lapply(names(vals_list), function(yr) {
      m <- vals_list[[yr]]
      data.table(year = as.integer(yr), month = as.integer(names(m)), value = unlist(m, use.names = FALSE))
    })
    dt2 <- rbindlist(rows)
    setnames(dt2, "value", varname)
    dt2
  }

  avg_temp <- make_clim("avg_temp", list(
    "2015"=c("6"=18.2,"8"=24.6,"9"=23.0,"10"=14.1), "2016"=c("6"=21.7,"8"=25.7,"9"=19.0,"10"=13.1),
    "2017"=c("6"=19.7,"10"=11.8), "2018"=c("6"=21.3,"8"=25.4,"9"=20.3,"10"=14.5),
    "2019"=c("6"=21.5,"8"=24.3,"9"=20.2,"10"=15.5), "2020"=c("6"=20.2,"8"=25.2,"9"=23.1,"10"=16.6),
    "2021"=c("6"=19.1,"8"=25.3,"9"=18.3,"10"=12.5)))
  min_temp <- make_clim("min_temp", list(
    "2015"=c("6"=12.6,"8"=16.7,"9"=14.9,"10"=8.6), "2016"=c("6"=13.0,"8"=17.2,"9"=11.4,"10"=6.4),
    "2017"=c("6"=12.7,"10"=5.7), "2018"=c("6"=13.6,"8"=17.2,"9"=12.3,"10"=8.2),
    "2019"=c("6"=14.7,"8"=15.9,"9"=12.1,"10"=8.0), "2020"=c("6"=12.8,"8"=16.0,"9"=14.8,"10"=9.9),
    "2021"=c("6"=12.0,"8"=16.3,"9"=11.4,"10"=5.7)))
  max_temp <- make_clim("max_temp", list(
    "2015"=c("6"=25.4,"8"=33.2,"9"=31.9,"10"=21.2), "2016"=c("6"=30.1,"8"=34.3,"9"=27.7,"10"=21.5),
    "2017"=c("6"=27.8,"10"=19.9), "2018"=c("6"=30.2,"8"=34.2,"9"=28.8,"10"=22.4),
    "2019"=c("6"=29.9,"8"=33.3,"9"=29.1,"10"=24.7), "2020"=c("6"=28.5,"8"=34.4,"9"=32.3,"10"=25.3),
    "2021"=c("6"=27.0,"8"=34.2,"9"=26.0,"10"=20.7)))
  precip <- make_clim("precipitation", list(
    "2015"=c("6"=152.6,"8"=21.2,"9"=2.7,"10"=21.7), "2016"=c("6"=16.1,"8"=4.8,"9"=32.3,"10"=2.4),
    "2017"=c("6"=91.2,"10"=33.4), "2018"=c("6"=22.3,"8"=52.6,"9"=2.0,"10"=67.2),
    "2019"=c("6"=89.6,"8"=9.3,"9"=1.0,"10"=13.8), "2020"=c("6"=62.5,"8"=0.4,"9"=8.9,"10"=39.2),
    "2021"=c("6"=118.8,"8"=4.4,"9"=10.4,"10"=25.5)))
  hot_days <- make_clim("hot_days", list(
    "2015"=c("6"=0,"8"=28,"9"=22,"10"=0), "2016"=c("6"=15,"8"=27,"9"=12,"10"=0),
    "2017"=c("6"=8,"10"=0), "2018"=c("6"=17,"8"=31,"9"=10,"10"=0),
    "2019"=c("6"=15,"8"=27,"9"=15,"10"=2), "2020"=c("6"=12,"8"=31,"9"=22,"10"=3),
    "2021"=c("6"=10,"8"=30,"9"=3,"10"=0)))
  avg_humidity <- make_clim("avg_humidity", list(
    "2015"=c("6"=73.8,"8"=47.3,"9"=43.4,"10"=68.4), "2016"=c("6"=53.4,"8"=43.3,"9"=50.4,"10"=56.7),
    "2017"=c("6"=63.2,"10"=61.7), "2018"=c("6"=59.6,"8"=43.9,"9"=48.5,"10"=68.7),
    "2019"=c("6"=69.9,"8"=49.9,"9"=48.9,"10"=60.4), "2020"=c("6"=68.5,"8"=42.8,"9"=51.1,"10"=64.0),
    "2021"=c("6"=71.0,"8"=47.7,"9"=61.3,"10"=67.3)))
  min_humidity <- make_clim("min_humidity", list(
    "2015"=c("6"=38.9,"8"=20.4,"9"=18.5,"10"=38.5), "2016"=c("6"=22.8,"8"=17.6,"9"=22.6,"10"=29.0),
    "2017"=c("6"=32.5,"10"=31.7), "2018"=c("6"=27.3,"8"=18.5,"9"=23.9,"10"=37.7),
    "2019"=c("6"=34.5,"8"=23.8,"9"=23.3,"10"=31.0), "2020"=c("6"=33.7,"8"=20.3,"9"=24.6,"10"=36.3),
    "2021"=c("6"=38.7,"8"=21.8,"9"=35.9,"10"=41.3)))
  max_humidity <- make_clim("max_humidity", list(
    "2015"=c("6"=98.9,"8"=77.2,"9"=71.2,"10"=91.9), "2016"=c("6"=89.1,"8"=73.2,"9"=80.6,"10"=82.8),
    "2017"=c("6"=89.0,"10"=85.3), "2018"=c("6"=89.8,"8"=72.1,"9"=75.9,"10"=91.4),
    "2019"=c("6"=96.6,"8"=78.0,"9"=76.4,"10"=88.1), "2020"=c("6"=94.9,"8"=69.9,"9"=80.6,"10"=86.4),
    "2021"=c("6"=97.7,"8"=76.0,"9"=87.2,"10"=91.4)))

  clim_all <- Reduce(function(a, b) merge(a, b, by = c("year","month"), all = TRUE),
                      list(avg_temp, min_temp, max_temp, precip, hot_days, avg_humidity, min_humidity, max_humidity))
  clim_all[, temp_amplitude := max_temp - min_temp]
  clim_all[, humidity_range := max_humidity - min_humidity]
  setorder(clim_all, month, year)
  clim_all[, delta_avg_temp  := avg_temp  - shift(avg_temp),  by = month]
  clim_all[, delta_min_temp  := min_temp  - shift(min_temp),  by = month]
  clim_all[, delta_max_temp  := max_temp  - shift(max_temp),  by = month]
  clim_all[, delta_humidity  := avg_humidity - shift(avg_humidity), by = month]
  clim_all[, delta_precip    := precipitation - shift(precipitation), by = month]

  CLIM_VARS <- c("avg_temp","min_temp","max_temp","precipitation","hot_days",
                 "avg_humidity","min_humidity","max_humidity","temp_amplitude","humidity_range",
                 "delta_avg_temp","delta_min_temp","delta_max_temp","delta_humidity","delta_precip")

  fst_meta <- merge(sample_meta[, .(sample_id, year, month)],
                     clim_all[, c("year","month",CLIM_VARS), with = FALSE],
                     by = c("year","month"), all.x = TRUE)
  setorder(fst_meta, year, month)

  clim_matrix <- as.matrix(fst_meta[, CLIM_VARS, with = FALSE])
  complete_rows <- complete.cases(clim_matrix)
  cat("Cohorts with complete climatic records (used for PCA fit + FST regression):",
      sum(complete_rows), "\n")

  pca_fit <- prcomp(clim_matrix[complete_rows, ], scale. = TRUE, center = TRUE)
  pct <- summary(pca_fit)$importance[2, ] * 100
  cat("PC1 + PC2 variance explained:", round(sum(pct[1:2]), 1), "% (PC1:",
      round(pct[1],1), "%, PC2:", round(pct[2],1), "%)\n")

  pc_scores <- as.data.table(pca_fit$x[, 1:2])
  setnames(pc_scores, c("PC1","PC2"))
  pc_scores[, sample_id := fst_meta$sample_id[complete_rows]]

  # --- Mean pairwise FST per cohort, from the 25-cohort (no 2015_09) sliding file ---
  pop_names_25 <- c("2015_06","2015_08","2015_10","2016_06","2016_08","2016_09","2016_10",
                     "2017_06","2017_10","2018_06","2018_08","2018_09","2018_10",
                     "2019_06","2019_08","2019_09","2019_10","2020_06","2020_08","2020_09","2020_10",
                     "2021_06","2021_08","2021_09","2021_10")

  fst_data <- fread(FST_FILE, data.table = FALSE, header = FALSE)
  sample_cols_idx <- 6:ncol(fst_data)
  pair_names <- sub("=.*", "", as.character(fst_data[1, sample_cols_idx]))
  fst_numeric <- apply(fst_data[, sample_cols_idx], 2, function(col) as.numeric(sub(".*=", "", as.character(col))))
  mean_fst <- colMeans(fst_numeric, na.rm = TRUE)
  names(mean_fst) <- pair_names

  fst_long <- data.table(i = as.integer(sub(":.*", "", names(mean_fst))),
                          j = as.integer(sub(".*:", "", names(mean_fst))),
                          FST = as.numeric(mean_fst))
  fst_mirror <- rbindlist(list(data.table(idx = fst_long$i, FST = fst_long$FST),
                                data.table(idx = fst_long$j, FST = fst_long$FST)))
  mean_pw_fst <- fst_mirror[, .(MeanFST = mean(FST, na.rm = TRUE)), by = idx]
  mean_pw_fst[, sample_id_25 := pop_names_25[idx]]

  # Map pop_names_25 (already "YYYY_MM" with zero-padded month, e.g. "2015_06")
  # directly onto sample_meta's own "YYYY_MM" label -- formats match exactly,
  # no string surgery needed.
  sample_meta[, ym_label := sprintf("%d_%02d", year, month)]
  setnames(mean_pw_fst, "sample_id_25", "ym_label")
  fst_by_sample <- merge(mean_pw_fst, sample_meta[, .(sample_id, ym_label)], by = "ym_label")

  climate_fst <- merge(pc_scores, fst_by_sample[, .(sample_id, MeanFST)], by = "sample_id")
  cat("Cohorts entering FST ~ PC1 + PC2 regression: n =", nrow(climate_fst), "\n")

  if (nrow(climate_fst) >= 10) {
    lm_pc <- lm(MeanFST ~ PC1 + PC2, data = climate_fst)
    sw_fst <- shapiro.test(residuals(lm_pc))

    sink(file.path(OUT_DIR, "fst_pc_assumption_summary.txt"))
    cat("FST ~ PC1 + PC2 model (n =", nrow(climate_fst), ")\n")
    print(summary(lm_pc))
    cat("\nShapiro-Wilk test on residuals:\n")
    print(sw_fst)
    sink()

    pdf(file.path(OUT_DIR, "fst_pc_residual_diagnostics.pdf"), width = 8, height = 6)
    par(mfrow = c(2, 2))
    plot(lm_pc)
    dev.off()

    cat("FST~PC1+PC2 Shapiro-Wilk p-value:", round(sw_fst$p.value, 4), "\n")
    cat(">>> Saved: fst_pc_assumption_summary.txt, fst_pc_residual_diagnostics.pdf\n")
  } else {
    cat("WARNING: fewer than 10 cohorts matched for FST~PC1+PC2 diagnostics; check ym_label merge.\n")
  }
}

# -----------------------------------------------------------------------
# 7. Diagnostic plot: pooled QQ-plot of residuals (random sample of SNPs)
# -----------------------------------------------------------------------
one_snp <- snps[snp_group == "random"][1]
sub_ex  <- y_master[chr == one_snp$chr & pos == one_snp$pos]
m_ex    <- lm(AF ~ year_val, data = sub_ex)

pdf(file.path(OUT_DIR, "example_residual_diagnostics.pdf"), width = 8, height = 6)
par(mfrow = c(2, 2))
plot(m_ex)
dev.off()

p_hist <- ggplot(diag_res, aes(x = shapiro_p, fill = snp_group)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0.05, linetype = "dashed", color = "red") +
  labs(title = "Distribution of Shapiro-Wilk p-values across sampled SNPs",
       x = "Shapiro-Wilk p-value (residuals of AF ~ Year)", y = "Count") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT_DIR, "shapiro_p_distribution.pdf"), p_hist, width = 7, height = 5)

cat("\n>>> All diagnostics done. Outputs in", OUT_DIR, "\n")
RSCRIPT

Rscript /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/ASSUMPTION_CHECK/assumption_check.R

echo ">>> Job finished: $(date)"
