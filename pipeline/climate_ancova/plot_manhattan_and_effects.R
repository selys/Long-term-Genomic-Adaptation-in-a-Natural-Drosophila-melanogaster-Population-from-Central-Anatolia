library(data.table)
library(ggplot2)
library(gridExtra)

BASE <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS"
OUT  <- BASE

# ============================================================
# 1. MANHATTAN PLOT
# ============================================================
cat(">>> Reading SNP data...\n")
sig   <- fread(file.path(BASE, "sig_year_only.tsv"))
full  <- fread(file.path(BASE, "snp_climate_ancova_full.tsv"))
hm    <- fread(file.path(BASE, "correct_high_moderate_unique.txt"),
               header=FALSE, col.names=c("chr","pos","gene","impact"))

# Chromosome order and cumulative positions
chr_order <- c("2L","2R","3L","3R","4")
chr_sizes  <- c("2L"=23513712, "2R"=25286936, "3L"=28110227, "3R"=32079331, "4"=1348131)

cum_offset <- c(0, cumsum(as.numeric(chr_sizes[-length(chr_sizes)])))
names(cum_offset) <- chr_order

full <- full[chr %in% chr_order]
full[, chr := factor(chr, levels=chr_order)]
full[, cum_pos := pos + cum_offset[as.character(chr)]]
full[, logp := -log10(p2_year)]
full[, sig := q2_year < 0.05]

# Merge HIGH/MODERATE flag
hm[, hm_flag := TRUE]
full <- merge(full, hm[, .(chr, pos, hm_flag)], by=c("chr","pos"), all.x=TRUE)
full[is.na(hm_flag), hm_flag := FALSE]

# Chromosome midpoints for x-axis labels
chr_mid <- cum_offset + chr_sizes/2

# Colors alternating per chromosome
chr_colors <- c("2L"="#4E79A7","2R"="#A0CBE8","3L"="#F28E2B","3R"="#FFBE7D","4"="#59A14F")

cat("Plotting Manhattan...\n")
p_man <- ggplot() +
  # Non-significant SNPs (grey, small)
  geom_point(data=full[sig==FALSE],
             aes(x=cum_pos, y=logp, color=chr),
             size=0.3, alpha=0.4) +
  # Year-significant SNPs
  geom_point(data=full[sig==TRUE & hm_flag==FALSE],
             aes(x=cum_pos, y=logp),
             color="#E15759", size=0.8, alpha=0.7) +
  # HIGH/MODERATE impact SNPs (highlighted)
  geom_point(data=full[hm_flag==TRUE],
             aes(x=cum_pos, y=logp),
             color="#000000", size=1.5, shape=18) +
  scale_color_manual(values=chr_colors, guide="none") +
  scale_x_continuous(breaks=chr_mid, labels=chr_order) +
  labs(
    title="Manhattan Plot — Year-Significant SNPs",
    subtitle="Yesiloz D. melanogaster (2015-2021) | Red: year-significant (q<0.05) | Black diamond: HIGH/MODERATE impact",
    x="Chromosome",
    y=expression(-log[10](p))
  ) +
  theme_minimal(base_size=12) +
  theme(
    panel.grid.minor=element_blank(),
    panel.grid.major.x=element_blank(),
    plot.subtitle=element_text(size=8, color="grey50")
  )

ggsave(file.path(OUT, "manhattan_year_sig.png"), p_man, width=14, height=5, dpi=150)
ggsave(file.path(OUT, "manhattan_year_sig.pdf"), p_man, width=14, height=5)
cat("Manhattan plot saved.\n")

# ============================================================
# 2. SNP EFFECT CATEGORY BAR PLOT
# ============================================================
cat(">>> Reading annotation data...\n")
ann <- fread(file.path(BASE, "year_sig_gene_annotations.txt"),
             header=TRUE,
             col.names=c("chr","pos","ref","alt","gene","effect","impact","feature"))

# Unique SNP x effect (collapse multi-transcript)
ann_u <- unique(ann[, .(chr, pos, effect, impact)])

# Simplify effect categories
ann_u[, effect_simple := fcase(
  grepl("missense", effect),                "Missense",
  grepl("synonymous", effect),              "Synonymous",
  grepl("stop_gained|stop_lost", effect),   "Stop gain/loss",
  grepl("splice", effect),                  "Splice region",
  grepl("5_prime_UTR|5prime", effect),      "5' UTR",
  grepl("3_prime_UTR|3prime", effect),      "3' UTR",
  grepl("intron", effect),                  "Intron",
  grepl("upstream", effect),               "Upstream",
  grepl("downstream", effect),             "Downstream",
  grepl("intergenic", effect),              "Intergenic",
  grepl("non_coding", effect),              "Non-coding exon",
  default = "Other"
)]

# Count unique SNPs per category (one SNP counted once per category)
eff_counts <- ann_u[, .N, by=effect_simple]
setorder(eff_counts, -N)

# Color by functional relevance
eff_counts[, category := fcase(
  effect_simple %in% c("Missense","Stop gain/loss","Splice region"), "High/Moderate",
  effect_simple %in% c("Synonymous","5' UTR","3' UTR"),              "Low",
  default = "Modifier"
)]

colors_eff <- c("High/Moderate"="#E15759", "Low"="#F28E2B", "Modifier"="#4E79A7")

p_eff <- ggplot(eff_counts, aes(x=reorder(effect_simple, N), y=N, fill=category)) +
  geom_bar(stat="identity", width=0.7) +
  geom_text(aes(label=N), hjust=-0.1, size=3.5) +
  scale_fill_manual(values=colors_eff, name="Functional impact") +
  coord_flip() +
  labs(
    title="Functional Effect Distribution of Year-Significant SNPs",
    subtitle="Yesiloz D. melanogaster (2015-2021) | n = 9,055 SNPs",
    x=NULL,
    y="Number of annotations"
  ) +
  theme_minimal(base_size=12) +
  theme(
    panel.grid.minor=element_blank(),
    legend.position="bottom",
    plot.subtitle=element_text(size=8, color="grey50")
  ) +
  expand_limits(y=max(eff_counts$N)*1.12)

ggsave(file.path(OUT, "snp_effect_distribution.png"), p_eff, width=10, height=7, dpi=150)
ggsave(file.path(OUT, "snp_effect_distribution.pdf"), p_eff, width=10, height=7)
cat("Effect distribution plot saved.\n")

# ============================================================
# 3. IMPACT SUMMARY (pie-like bar)
# ============================================================
imp_counts <- ann_u[, .N, by=impact]
imp_counts <- imp_counts[impact %in% c("HIGH","MODERATE","LOW","MODIFIER")]
setorder(imp_counts, -N)

imp_colors <- c("HIGH"="#B22222","MODERATE"="#E15759","LOW"="#F28E2B","MODIFIER"="#4E79A7")

p_imp <- ggplot(imp_counts, aes(x=reorder(impact, N), y=N, fill=impact)) +
  geom_bar(stat="identity", width=0.6) +
  geom_text(aes(label=format(N, big.mark=",")), hjust=-0.1, size=4) +
  scale_fill_manual(values=imp_colors, guide="none") +
  coord_flip() +
  labs(
    title="Functional Impact Distribution of Year-Significant SNPs",
    subtitle="Yesiloz D. melanogaster (2015-2021) | n = 9,055 SNPs",
    x=NULL, y="Number of annotations"
  ) +
  theme_minimal(base_size=12) +
  theme(panel.grid.minor=element_blank(),
        plot.subtitle=element_text(size=8, color="grey50")) +
  expand_limits(y=max(imp_counts$N)*1.12)

ggsave(file.path(OUT, "snp_impact_distribution.png"), p_imp, width=8, height=4, dpi=150)
ggsave(file.path(OUT, "snp_impact_distribution.pdf"), p_imp, width=8, height=4)
cat("Impact distribution plot saved.\n")

cat("\n>>> All plots saved to", OUT, "\n")
