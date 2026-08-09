#!/bin/bash
#SBATCH -p hamsi
#SBATCH -A ssenkal
#SBATCH -J fst_time_mantel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=8G
#SBATCH -t 01:00:00
#SBATCH --qos=normal
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/fst_time_mantel_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/fst_time_mantel_%j.err

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

OUT_DIR="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/FST_TIME_DEPENDENCE"
mkdir -p $OUT_DIR

FST_FILE="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/fst_results/mega_all_years_25_sliding_fst_no2015_09.txt"
if [ ! -f "$FST_FILE" ]; then
    echo "HATA: $FST_FILE bulunamadi."
    exit 1
fi
echo "Using: $FST_FILE"

cat > $OUT_DIR/fst_time_mantel.R << 'RSCRIPT'
# =============================================================================
# Is pairwise FST time-dependent? (Reviewer 2: "this is not visible from
# figure 1 and there is no statistical test")
#
# Approach: Mantel test between the pairwise FST matrix (25 cohorts,
# 2015_09 excluded) and a pairwise TEMPORAL DISTANCE matrix (elapsed time
# between cohorts, in years). Implemented as a permutation test from
# scratch (no 'vegan'/'ade4' dependency) -- rows/columns of the temporal
# distance matrix are randomly relabeled many times to build a null
# distribution for the observed FST~time-distance correlation.
# =============================================================================
library(data.table)
library(ggplot2)

set.seed(42)
OUT_DIR <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/FST_TIME_DEPENDENCE"
FST_FILE <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/fst_results/mega_all_years_25_sliding_fst_no2015_09.txt"

# -----------------------------------------------------------------------
# 1. Cohort order (25 cohorts, 2015_09 excluded) -- matches the order the
#    sliding-window FST file's samples were provided in.
# -----------------------------------------------------------------------
pop_names_25 <- c("2015_06","2015_08","2015_10","2016_06","2016_08","2016_09","2016_10",
                   "2017_06","2017_10","2018_06","2018_08","2018_09","2018_10",
                   "2019_06","2019_08","2019_09","2019_10","2020_06","2020_08","2020_09","2020_10",
                   "2021_06","2021_08","2021_09","2021_10")
n_pops <- length(pop_names_25)

year_of  <- as.integer(substr(pop_names_25, 1, 4))
month_of <- as.integer(substr(pop_names_25, 6, 7))
decimal_time <- year_of + (month_of - 1) / 12   # continuous time axis, in years

# -----------------------------------------------------------------------
# 2. Build the pairwise FST matrix (same parsing logic as the heatmap /
#    PCA scripts used throughout this project)
# -----------------------------------------------------------------------
fst_data <- fread(FST_FILE, data.table = FALSE, header = FALSE)
sample_cols_idx <- 6:ncol(fst_data)
pair_names <- sub("=.*", "", as.character(fst_data[1, sample_cols_idx]))
fst_numeric <- apply(fst_data[, sample_cols_idx], 2, function(col) as.numeric(sub(".*=", "", as.character(col))))
mean_fst <- colMeans(fst_numeric, na.rm = TRUE)
names(mean_fst) <- pair_names

i_idx <- as.integer(sub(":.*", "", names(mean_fst)))
j_idx <- as.integer(sub(".*:", "", names(mean_fst)))

FST <- matrix(NA_real_, n_pops, n_pops)
for (k in seq_along(mean_fst)) {
  FST[i_idx[k], j_idx[k]] <- mean_fst[k]
  FST[j_idx[k], i_idx[k]] <- mean_fst[k]
}
diag(FST) <- 0
cat("FST matrix built:", n_pops, "x", n_pops, "cohorts\n")
cat("Missing (NA) pairwise entries:", sum(is.na(FST[upper.tri(FST)])), "of", sum(upper.tri(FST)), "\n")

# -----------------------------------------------------------------------
# 3. Build the temporal distance matrix (absolute time elapsed, in years)
# -----------------------------------------------------------------------
TIME_DIST <- as.matrix(dist(decimal_time, method = "euclidean"))
rownames(TIME_DIST) <- colnames(TIME_DIST) <- pop_names_25
rownames(FST) <- colnames(FST) <- pop_names_25

# -----------------------------------------------------------------------
# 4. Mantel test (permutation-based)
# -----------------------------------------------------------------------
mantel_test <- function(mat1, mat2, n_perm = 9999) {
  ut <- upper.tri(mat1)
  valid <- ut & !is.na(mat1) & !is.na(mat2)
  v1 <- mat1[valid]
  v2 <- mat2[valid]
  obs_r <- cor(v1, v2, method = "pearson")

  n <- nrow(mat1)
  perm_r <- numeric(n_perm)
  for (p in seq_len(n_perm)) {
    perm_order <- sample(n)
    mat2_perm <- mat2[perm_order, perm_order]
    v2_perm <- mat2_perm[valid]
    perm_r[p] <- cor(v1, v2_perm, method = "pearson")
  }
  # one-sided p-value: is the observed correlation greater than expected by chance?
  p_value <- (sum(perm_r >= obs_r) + 1) / (n_perm + 1)
  list(obs_r = obs_r, perm_r = perm_r, p_value = p_value, n_pairs = sum(valid))
}

cat("\n>>> Running Mantel test (9,999 permutations)...\n")
mt <- mantel_test(FST, TIME_DIST, n_perm = 9999)
cat(sprintf("Mantel r = %.4f, p = %.4g (one-sided, n pairs = %d)\n",
            mt$obs_r, mt$p_value, mt$n_pairs))

# -----------------------------------------------------------------------
# 5. Complementary descriptive regression (effect size / R^2 for reporting;
#    the Mantel permutation p-value above is the valid significance test --
#    this lm is for interpretability only, since pairwise FST values are
#    not independent observations and an ordinary lm p-value would be
#    invalid/anti-conservative here).
# -----------------------------------------------------------------------
ut <- upper.tri(FST)
valid <- ut & !is.na(FST)
df_pairs <- data.frame(
  FST = FST[valid],
  time_dist = TIME_DIST[valid]
)
lm_fit <- lm(FST ~ time_dist, data = df_pairs)
lm_r2 <- summary(lm_fit)$r.squared
lm_slope <- coef(lm_fit)[["time_dist"]]

# -----------------------------------------------------------------------
# 6. Save results
# -----------------------------------------------------------------------
sink(file.path(OUT_DIR, "fst_time_mantel_summary.txt"))
cat("Is pairwise FST time-dependent? Mantel test\n")
cat("=============================================\n\n")
cat("Cohorts:", n_pops, "(2015_09 excluded)\n")
cat("Valid pairwise comparisons:", mt$n_pairs, "\n\n")
cat("Mantel test (FST distance matrix vs. temporal distance matrix):\n")
cat(sprintf("  Observed Mantel r = %.4f\n", mt$obs_r))
cat(sprintf("  Permutation p-value (n=9,999, one-sided) = %.4g\n\n", mt$p_value))
cat("Descriptive regression (FST ~ temporal distance in years) -- for effect\n")
cat("size / interpretability only; the Mantel p-value above is the valid\n")
cat("significance test given non-independence among pairwise FST values:\n")
cat(sprintf("  slope = %.5f FST per year of separation\n", lm_slope))
cat(sprintf("  R^2 = %.4f\n", lm_r2))
sink()

fwrite(df_pairs, file.path(OUT_DIR, "fst_vs_temporal_distance_pairs.tsv"), sep = "\t")

p1 <- ggplot(df_pairs, aes(x = time_dist, y = FST)) +
  geom_point(alpha = 0.5, size = 1.8) +
  geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
  labs(title = "Pairwise FST increases with temporal distance between cohorts",
       subtitle = sprintf("Mantel r = %.3f, permutation p = %.3g (n = %d pairs, 9,999 permutations)",
                           mt$obs_r, mt$p_value, mt$n_pairs),
       x = "Temporal distance between cohorts (years)",
       y = expression(Pairwise~italic(F)[ST])) +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT_DIR, "fst_vs_temporal_distance.pdf"), p1, width = 7.5, height = 5.5)
ggsave(file.path(OUT_DIR, "fst_vs_temporal_distance.png"), p1, width = 7.5, height = 5.5, dpi = 300)

cat("\n>>> Done. Outputs in", OUT_DIR, "\n")
RSCRIPT

Rscript $OUT_DIR/fst_time_mantel.R

echo ">>> Job finished: $(date)"
