library(data.table)
library(ggplot2)

OUT_DIR <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL"

# --- Gene lengths from GTF/GFF3 ---
# NOTE: the reference genome here is NCBI RefSeq (GCF_000001215.4), so the
# attribute column may be GFF3-style ("gene=Ank2;...", "Name=Ank2;...")
# rather than GTF-style ('gene_name "Ank2";'). Try both, in order of
# preference, instead of assuming GTF quoting.
gtf <- fread(file.path(OUT_DIR, "genes_raw.gtf"), header = FALSE, sep = "\t", quote = "")
setnames(gtf, c("chr","source","feature","start","end","score","strand","frame","attr"))
gtf[, length_bp := end - start + 1]

extract_attr <- function(attr, key) {
  gtf_pat <- paste0(key, ' "([^"]+)"')          # GTF:  key "value";
  gff_pat <- paste0('(^|;)\\s*', key, '=([^;]+)') # GFF3: key=value;
  out <- rep(NA_character_, length(attr))
  hit_gtf <- grepl(gtf_pat, attr)
  out[hit_gtf] <- sub(paste0('.*', gtf_pat, '.*'), '\\1', attr[hit_gtf])
  hit_gff <- is.na(out) & grepl(gff_pat, attr)
  out[hit_gff] <- sub(paste0('.*', gff_pat, '.*'), '\\2', attr[hit_gff])
  out
}

# Try, in order: gene_name (GTF) -> gene (GFF3, most common on RefSeq) ->
# Name (GFF3 fallback) -> gene_id / ID (last resort, may be FBgn/LOC IDs)
gtf[, gene_name := extract_attr(attr, "gene_name")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "gene")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "Name")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "gene_id")]
gtf[is.na(gene_name), gene_name := extract_attr(attr[is.na(gene_name)], "ID")]

n_unparsed <- sum(is.na(gtf$gene_name))
cat("Gene records with NO parsable name (dropped):", n_unparsed, "of", nrow(gtf), "\n")
if (n_unparsed > 0) {
  cat("Example unparsed attribute strings:\n")
  print(head(gtf[is.na(gene_name), attr], 3))
}
gtf <- gtf[!is.na(gene_name)]

# --- Restrict to autosomes only ---
# Uses the SAME 5 RefSeq accessions the main pipeline's mpileup filtering
# step used (2L=NT_033779.5, 2R=NT_033778.4, 3L=NT_037436.4, 3R=NT_033777.3,
# 4=NC_004353.4). X/Y/mitochondrial genes are excluded so this analysis
# covers exactly the same SNP universe as the rest of the manuscript
# ("SNPs on all autosomes").
AUTOSOME_ACCESSIONS <- c("NT_033779.5","NT_033778.4","NT_037436.4","NT_033777.3","NC_004353.4")
n_before <- nrow(gtf)
gtf <- gtf[chr %in% AUTOSOME_ACCESSIONS]
cat("Gene records restricted to autosomes:", nrow(gtf), "of", n_before,
    "(dropped X/Y/mito/unplaced scaffolds)\n")
if (nrow(gtf) == 0) {
  cat("WARNING: 0 rows after autosome filter -- the GFF's chromosome/seqid naming\n")
  cat("does not match the expected accessions. Example chr values seen in GTF:\n")
  print(unique(fread(file.path(OUT_DIR, "genes_raw.gtf"), header = FALSE, sep = "\t", quote = "")$V1))
}

gene_lengths <- gtf[, .(length_bp = max(length_bp)), by = gene_name]  # collapse duplicate records per gene
cat("Unique autosomal genes with length info:", nrow(gene_lengths), "\n")

# --- Background & significant SNP-gene tables ---
bg  <- fread(file.path(OUT_DIR, "background_chr_pos_gene.txt"), header = FALSE,
             col.names = c("chr","pos","gene_name"))
sig <- fread(file.path(OUT_DIR, "sig_chr_pos_gene.txt"), header = FALSE,
             col.names = c("chr","pos","gene_name"))

if (nrow(bg) == 0) {
  stop("background_chr_pos_gene.txt loaded with 0 rows -- the bcftools extraction ",
       "step failed upstream (see the .sh script's own check, which should have ",
       "already caught this before R even started). Re-run the bash script and ",
       "inspect its output/error log.")
}

# --- Sanity check BEFORE merging: do gene symbols actually overlap? ---
# This is the check that would have caught a gene_name/gene_id mismatch
# (e.g. VCF giving "Ank2" while the GTF only yielded "FBgn0000083") before
# silently producing a near-empty merged table.
vcf_genes <- unique(bg$gene_name)
overlap_n <- length(intersect(vcf_genes, gene_lengths$gene_name))
overlap_pct <- round(100 * overlap_n / length(vcf_genes), 1)
cat(sprintf("\nGene symbol overlap check: %d / %d (%.1f%%) of VCF-annotated genes found in GTF length table\n",
            overlap_n, length(vcf_genes), overlap_pct))
if (overlap_pct < 50) {
  cat("WARNING: overlap below 50% -- gene_name/gene_id mismatch is likely.\n")
  cat("Example VCF gene symbols NOT found in GTF table:\n")
  print(head(setdiff(vcf_genes, gene_lengths$gene_name), 10))
  cat("Example GTF gene_name values (for comparison):\n")
  print(head(gene_lengths$gene_name, 10))
  cat("Fix extract_attr()/attribute key choice above before trusting downstream results.\n")
}

bg_counts  <- bg[,  .(n_tested = .N), by = gene_name]
sig_counts <- sig[, .(n_sig    = .N), by = gene_name]

merged <- merge(bg_counts, sig_counts, by = "gene_name", all.x = TRUE)
merged[is.na(n_sig), n_sig := 0]
merged <- merge(merged, gene_lengths, by = "gene_name")
merged <- merged[length_bp > 0 & n_tested > 0]
merged[, prop_sig := n_sig / n_tested]

cat("\nGenes in final merged table:", nrow(merged), "\n")

# --- Test 1: raw sig-SNP count vs gene length (expected to correlate trivially) ---
ct1 <- cor.test(merged$length_bp, merged$n_sig, method = "spearman")
cat("\n== Raw significant-SNP count vs gene length (Spearman) ==\n")
print(ct1)

# --- Test 2: proportion significant (normalized by tested SNPs) vs gene length ---
#     This is the key test: if length no longer predicts the PROPORTION,
#     the raw-count correlation is just a density artifact, not a real bias.
ct2 <- cor.test(merged$length_bp, merged$prop_sig, method = "spearman")
cat("\n== Proportion significant (n_sig/n_tested) vs gene length (Spearman) ==\n")
print(ct2)

# --- Test 3: Poisson regression, n_sig ~ log(length), offset = log(n_tested) ---
#     Directly asks: controlling for how many SNPs were even testable in the
#     gene (which itself scales with length), does length still matter?
merged[, log_length := log(length_bp)]
pois_fit <- glm(n_sig ~ log_length + offset(log(n_tested)), data = merged, family = poisson())
cat("\n== Poisson regression: n_sig ~ log(gene length), offset = log(n_tested) ==\n")
print(summary(pois_fit))

# --- Save table and plots ---
fwrite(merged, file.path(OUT_DIR, "gene_length_vs_significance.tsv"), sep = "\t")

p1 <- ggplot(merged, aes(x = length_bp, y = n_sig)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_x_log10() +
  labs(title = "Raw significant SNP count vs gene length",
       x = "Gene length (bp, log10 scale)", y = "Year-significant SNPs per gene") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT_DIR, "raw_count_vs_length.pdf"), p1, width = 7, height = 5)

p2 <- ggplot(merged, aes(x = length_bp, y = prop_sig)) +
  geom_point(alpha = 0.3, size = 1) +
  scale_x_log10() +
  labs(title = "Proportion significant (normalized) vs gene length",
       x = "Gene length (bp, log10 scale)", y = "Proportion of tested SNPs that are year-significant") +
  theme_minimal(base_size = 12)
ggsave(file.path(OUT_DIR, "proportion_vs_length.pdf"), p2, width = 7, height = 5)

sink(file.path(OUT_DIR, "gene_length_control_summary.txt"))
cat("Gene-length confound check\n===========================\n\n")
cat("Genes analyzed:", nrow(merged), "\n\n")
cat("Test 1 - raw count vs length (Spearman rho, p):\n")
cat(sprintf("  rho = %.4f, p = %.3g\n\n", ct1$estimate, ct1$p.value))
cat("Test 2 - proportion significant vs length (Spearman rho, p):\n")
cat(sprintf("  rho = %.4f, p = %.3g\n\n", ct2$estimate, ct2$p.value))
cat("Test 3 - Poisson regression coefficient on log(length), offset=log(n_tested):\n")
print(coef(summary(pois_fit)))
sink()

cat("\n>>> Done. Outputs in", OUT_DIR, "\n")
