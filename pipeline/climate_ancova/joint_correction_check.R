#!/usr/bin/env Rscript
# joint_correction_check.R
# Verifies whether pooling BH correction across all Model 1-4 tests
# changes the year-significant SNP count reported in the manuscript.

library(data.table)

IN_FILE  <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/snp_climate_ancova_full.tsv"
OUT_DIR  <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS"

cat(">>> Reading:", IN_FILE, "\n")
res <- fread(IN_FILE)

# ---- 1. Pool every p-value column from Models 1-4 into one long vector ----
p_cols <- c("p1_season", "p2_year",
            "p3_season", "p3_year", "p3_interaction",
            "p4_season", "p4_year", "p4_PC1", "p4_PC2")

pooled <- rbindlist(lapply(p_cols, function(cn) {
  data.table(chr = res$chr, pos = res$pos, test = cn, pval = res[[cn]])
}))
pooled <- pooled[!is.na(pval)]
cat("Total pooled tests:", nrow(pooled), "\n")

# ---- 2. Apply BH jointly across ALL tests from ALL models ----
pooled[, q_joint := p.adjust(pval, method = "BH")]

# ---- 3. Original (separate, per-model) significant count for Model 2 (year) ----
n_orig_year <- res[q2_year < 0.05, .N]

# ---- 4. Joint-corrected significant count for Model 2 (year) ----
year_joint <- pooled[test == "p2_year" & q_joint < 0.05]
n_joint_year <- nrow(year_joint)

# ---- 5. Same check for season (Model 1) as a sanity comparison ----
n_orig_season   <- res[q1_season < 0.05, .N]
season_joint    <- pooled[test == "p1_season" & q_joint < 0.05]
n_joint_season  <- nrow(season_joint)

cat("\n============ RESULTS ============\n")
cat(sprintf("Year-significant SNPs   (separate correction, original): %d\n", n_orig_year))
cat(sprintf("Year-significant SNPs   (joint correction, all models pooled): %d\n", n_joint_year))
cat(sprintf("Season-significant SNPs (separate correction, original): %d\n", n_orig_season))
cat(sprintf("Season-significant SNPs (joint correction, all models pooled): %d\n", n_joint_season))

fwrite(year_joint[, .(chr, pos, pval, q_joint)],
       file.path(OUT_DIR, "year_significant_SNPs_joint_correction.tsv"), sep = "\t")

cat("\n>>> Saved joint-correction year-significant SNP list.\n")
cat(">>> Done:", format(Sys.time()), "\n")
