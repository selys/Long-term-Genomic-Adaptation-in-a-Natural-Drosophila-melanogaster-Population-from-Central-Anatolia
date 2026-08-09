#!/usr/bin/env Rscript
# =============================================================================
# Candidate gene density check
# Question: are the manuscript's "top candidate genes" (Ank2, Muc68Ca, etc.)
# genuinely enriched for year-significant SNPs, or do they just carry a lot
# of SNPs because they are large? Ranks each candidate gene's SNP DENSITY
# (not raw count) against the density distribution of all ~14,000 genes.
#
# Run directly (no sbatch needed, input table is small):
#   Rscript candidate_gene_density_check.R
# =============================================================================
library(data.table)

IN_FILE <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL/gene_length_vs_significance.tsv"
OUT_DIR <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL"

dt <- fread(IN_FILE)
cat("Genes loaded:", nrow(dt), "\n")

# Two density measures:
#  - prop_sig    = n_sig / n_tested   (already in the file: controls for how
#                  many SNPs were even callable in that gene)
#  - density_bp  = n_sig / length_bp  (SNPs per bp of gene length -- the
#                  literal "per-length" density the reviewer's comment and
#                  our discussion concerns)
dt[, density_bp := n_sig / length_bp]

dt[, pct_prop_sig   := rank(prop_sig,   ties.method = "average") / .N * 100]
dt[, pct_density_bp := rank(density_bp, ties.method = "average") / .N * 100]
dt[, pct_length     := rank(length_bp,  ties.method = "average") / .N * 100]

candidates <- c("Ank2", "Muc68Ca", "Muc68E", "trp", "dpy", "Mur29B", "Msp300")

result <- dt[gene_name %in% candidates]
missing <- setdiff(candidates, dt$gene_name)
if (length(missing) > 0) {
  cat("WARNING: not found in table (check spelling/case):", paste(missing, collapse=", "), "\n")
}

setorder(result, -pct_density_bp)

cat("\n================ Candidate gene enrichment check ================\n")
cat("pct_length     = size percentile among all genes (100 = largest)\n")
cat("pct_prop_sig   = percentile of (n_sig/n_tested), controls for testability\n")
cat("pct_density_bp = percentile of (n_sig/length_bp), raw per-bp density\n")
cat("If pct_density_bp/pct_prop_sig is near pct_length -> gene is 'big, not special'\n")
cat("If pct_density_bp/pct_prop_sig is well ABOVE pct_length -> genuine enrichment\n\n")

print(result[, .(gene_name, length_bp, n_tested, n_sig, prop_sig,
                  pct_length = round(pct_length,1),
                  pct_prop_sig = round(pct_prop_sig,1),
                  pct_density_bp = round(pct_density_bp,1))])

fwrite(result, file.path(OUT_DIR, "candidate_gene_density_check.tsv"), sep = "\t")

cat("\n>>> Saved: candidate_gene_density_check.tsv\n")
