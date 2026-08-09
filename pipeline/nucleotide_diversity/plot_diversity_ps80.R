library(data.table)
library(ggplot2)
library(gridExtra)

BASE <- "/arf/scratch/ssenkal/selin/sekanslar"
DIR  <- file.path(BASE, "all_years_26/PI_ANALYSIS_1k_ps80")

samples_list <- c(
  "2015_06","2015_08","2015_10",
  "2016_06","2016_08","2016_09","2016_10",
  "2017_06","2017_10",
  "2018_06","2018_08","2018_09","2018_10",
  "2019_06","2019_08","2019_09","2019_10",
  "2020_06","2020_08","2020_09","2020_10",
  "2021_06","2021_08","2021_09","2021_10"
)
# NOTE: 2015_09 intentionally excluded (low coverage, <8% valid windows)

read_div <- function(measure) {
  rbindlist(lapply(samples_list, function(s) {
    f <- file.path(DIR, paste0(s, "_", measure, ".txt"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f, header=FALSE, col.names=c("chr","pos","count","frac","mean"))
    dt <- dt[mean != "na" & !is.na(suppressWarnings(as.numeric(mean)))]
    dt[, mean := as.numeric(mean)]
    dt <- dt[frac >= 0.5]
    if (nrow(dt) == 0) return(NULL)
    data.table(sample=s, mean_val=mean(dt$mean, na.rm=TRUE), n_windows=nrow(dt))
  }), fill=TRUE)
}


read_div_D <- function(measure_suffix, dir) {
  rbindlist(lapply(samples_list, function(s) {
    f <- file.path(dir, paste0(s, "_", measure_suffix, ".txt"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f, header=FALSE, col.names=c("chr","pos","count","frac","mean"))
    dt <- dt[mean != "na" & !is.na(suppressWarnings(as.numeric(mean)))]
    dt[, mean := as.numeric(mean)]
    dt <- dt[frac >= 0.5]
    if (nrow(dt) == 0) return(NULL)
    data.table(sample=s, mean_val=mean(dt$mean, na.rm=TRUE), n_windows=nrow(dt))
  }), fill=TRUE)
}

pi_dat    <- read_div("pi")
theta_dat <- read_div("theta")
D_dat     <- read_div_D("D_ps80", file.path(DIR, "verify_ps80"))


cat("Pi windows per sample:\n"); print(pi_dat[, .(sample, n_windows)])

# Add month/year
add_meta <- function(dt) {
  dt[, year  := as.integer(substr(sample, 1, 4))]
  dt[, month := substr(sample, 6, 7)]
  dt[, month_label := factor(month,
      levels=c("06","08","09","10"),
      labels=c("June","August","September","October"))]
  dt
}
pi_dat    <- add_meta(pi_dat)
theta_dat <- add_meta(theta_dat)
D_dat     <- add_meta(D_dat)

# Save summaries
fwrite(pi_dat[,.(sample,year,month,mean_pi=mean_val,n_windows)],
       file.path(DIR,"pi_summary.tsv"), sep="\t")
fwrite(theta_dat[,.(sample,year,month,mean_theta=mean_val,n_windows)],
       file.path(DIR,"theta_summary.tsv"), sep="\t")
fwrite(D_dat[,.(sample,year,month,mean_D=mean_val,n_windows)],
       file.path(DIR,"D_summary.tsv"), sep="\t")

colors <- c("June"="#2ecc71","August"="#e74c3c",
            "September"="#f39c12","October"="#8e44ad")
sub <- "Yesiloz D. melanogaster (2015-2021) - 1kb windows, pool-size 80, autosomal only\n(2015_09 excluded: low coverage)"

make_plot <- function(dt, y_col, y_lab, title_str) {
  ggplot(dt, aes(x=sample, y=mean_val, group=1)) +
    geom_line(color="steelblue", linewidth=0.7) +
    geom_point(aes(color=month_label), size=3.5) +
    scale_color_manual(values=colors, name="Month") +
    scale_x_discrete(guide=guide_axis(angle=45)) +
    labs(title=title_str, subtitle=sub, x="Sample", y=y_lab) +
    theme_minimal(base_size=11) +
    theme(plot.subtitle=element_text(size=8, color="grey50"),
          axis.text.x=element_text(size=7))
}

p_pi    <- make_plot(pi_dat,    "mean_val", "Mean pi",    "Nucleotide Diversity (pi)")
p_theta <- make_plot(theta_dat, "mean_val", "Mean theta", "Watterson's Theta (theta)")
p_D     <- make_plot(D_dat,     "mean_val", "Mean D",     "Tajima's D")

# Individual PNGs
ggsave(file.path(DIR,"pi_temporal.png"),    p_pi,    width=14, height=5, dpi=150)
ggsave(file.path(DIR,"theta_temporal.png"), p_theta, width=14, height=5, dpi=150)
ggsave(file.path(DIR,"D_temporal.png"),     p_D,     width=14, height=5, dpi=150)

# Individual PDFs
ggsave(file.path(DIR,"pi_temporal.pdf"),    p_pi,    width=14, height=5)
ggsave(file.path(DIR,"theta_temporal.pdf"), p_theta, width=14, height=5)
ggsave(file.path(DIR,"D_temporal.pdf"),     p_D,     width=14, height=5)

# Combined PDF
pdf(file.path(DIR,"diversity_all_temporal.pdf"), width=14, height=13)
grid.arrange(p_pi, p_theta, p_D, ncol=1)
dev.off()

cat(">>> Done! Plots saved to", DIR, "\n")
