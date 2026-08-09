library(data.table)
library(ggplot2)

BASE <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS"
OUT  <- BASE

cat(">>> Reading annotation data...\n")
ann <- fread(file.path(BASE, "year_sig_gene_annotations.txt"),
             header=TRUE,
             col.names=c("chr","pos","ref","alt","gene","effect","impact","feature"))

cat("Total annotation rows:", nrow(ann), "\n")

# ============================================================
# Her SNP (chr+pos) için unique effect kategorilerini al
# Aynı SNP'in aynı effect'i birden fazla transcriptte varsa bir kez say
# ============================================================
ann_u <- unique(ann[, .(chr, pos, effect, impact)])
cat("Unique chr+pos+effect+impact combinations:", nrow(ann_u), "\n")
cat("Unique SNPs:", uniqueN(ann_u[, .(chr, pos)]), "\n")

# Effect kategorilerini sadeleştir
ann_u[, effect_simple := fcase(
  grepl("missense", effect),              "Missense",
  grepl("stop_gained|stop_lost", effect), "Stop gain/loss",
  grepl("splice", effect),               "Splice region",
  grepl("synonymous", effect),           "Synonymous",
  grepl("5_prime_UTR|5prime", effect),   "5' UTR",
  grepl("3_prime_UTR|3prime", effect),   "3' UTR",
  grepl("intron", effect),               "Intron",
  grepl("upstream", effect),             "Upstream",
  grepl("downstream", effect),           "Downstream",
  grepl("intergenic", effect),           "Intergenic",
  grepl("non_coding", effect),           "Non-coding exon",
  default = "Other"
)]

# Her SNP için unique effect_simple kategorilerini say
# (bir SNP aynı kategoride birden fazla gene için varsa bir kez say)
ann_snp_eff <- unique(ann_u[, .(chr, pos, effect_simple, impact)])

# Kategori bazında SNP sayısı
eff_counts <- ann_snp_eff[, .(n_snps = uniqueN(paste(chr, pos))), by=effect_simple]
setorder(eff_counts, -n_snps)

cat("\nEffect distribution (unique SNPs per category):\n")
print(eff_counts)
cat("Note: total > 9055 because some SNPs fall in multiple categories\n")

# Renk — impact grubuna göre
eff_counts[, category := fcase(
  effect_simple %in% c("Missense","Stop gain/loss","Splice region"), "High/Moderate",
  effect_simple %in% c("Synonymous","5' UTR","3' UTR"),              "Low",
  default = "Modifier"
)]

colors_eff <- c("High/Moderate"="#E15759", "Low"="#F28E2B", "Modifier"="#4E79A7")

p_eff <- ggplot(eff_counts, aes(x=reorder(effect_simple, n_snps), y=n_snps, fill=category)) +
  geom_bar(stat="identity", width=0.7) +
  geom_text(aes(label=n_snps), hjust=-0.1, size=3.5) +
  scale_fill_manual(values=colors_eff, name="Functional impact") +
  coord_flip() +
  labs(
    title="Functional Effect Distribution of Year-Significant SNPs",
    subtitle="Yesiloz D. melanogaster (2015-2021) | n = 9,055 SNPs\n(bars sum to >9,055 as some SNPs have multiple annotations across genes)",
    x=NULL,
    y="Number of SNPs"
  ) +
  theme_minimal(base_size=12) +
  theme(
    panel.grid.minor=element_blank(),
    legend.position="bottom",
    plot.subtitle=element_text(size=8, color="grey50")
  ) +
  expand_limits(y=max(eff_counts$n_snps)*1.15)

ggsave(file.path(OUT, "snp_effect_distribution_v2.png"), p_eff, width=10, height=7, dpi=150)
ggsave(file.path(OUT, "snp_effect_distribution_v2.pdf"), p_eff, width=10, height=7)
cat("Effect distribution plot saved.\n")

# ============================================================
# IMPACT dağılımı da aynı mantıkla
# ============================================================
imp_snp <- unique(ann_u[, .(chr, pos, impact)])
imp_counts <- imp_snp[impact %in% c("HIGH","MODERATE","LOW","MODIFIER"),
                      .(n_snps = uniqueN(paste(chr, pos))), by=impact]
setorder(imp_counts, -n_snps)

cat("\nImpact distribution (unique SNPs per category):\n")
print(imp_counts)

imp_colors <- c("HIGH"="#B22222","MODERATE"="#E15759","LOW"="#F28E2B","MODIFIER"="#4E79A7")

p_imp <- ggplot(imp_counts, aes(x=reorder(impact, n_snps), y=n_snps, fill=impact)) +
  geom_bar(stat="identity", width=0.6) +
  geom_text(aes(label=format(n_snps, big.mark=",")), hjust=-0.1, size=4) +
  scale_fill_manual(values=imp_colors, guide="none") +
  coord_flip() +
  labs(
    title="Functional Impact Distribution of Year-Significant SNPs",
    subtitle="Yesiloz D. melanogaster (2015-2021) | n = 9,055 SNPs\n(bars sum to >9,055 as some SNPs have multiple annotations across genes)",
    x=NULL, y="Number of SNPs"
  ) +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(),
        plot.subtitle=element_text(size=8, color="grey50")) +
  expand_limits(y=max(imp_counts$n_snps)*1.15)

ggsave(file.path(OUT, "snp_impact_distribution_v2.png"), p_imp, width=8, height=4, dpi=150)
ggsave(file.path(OUT, "snp_impact_distribution_v2.pdf"), p_imp, width=8, height=4)
cat("Impact distribution plot saved.\n")

cat("\n>>> All done!\n")
