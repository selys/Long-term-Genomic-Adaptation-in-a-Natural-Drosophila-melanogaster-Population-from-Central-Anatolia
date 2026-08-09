library(data.table)
library(ggplot2)
library(gridExtra)
library(readxl)

BASE <- "/arf/scratch/ssenkal/selin/sekanslar"
DIR  <- file.path(BASE, "all_years_26/PI_ANALYSIS_1k_ps80")

# ============================================================
# 1. WINTER SEVERITY SCORES
# ============================================================
don     <- fread(file.path(BASE, "don_gunler.csv"))
soguk   <- fread(file.path(BASE, "soguk_gunler.csv"))
mintemp <- fread(file.path(BASE, "min_sicaklik.csv"))

kis_skoru <- rbindlist(lapply(2015:2021, function(y) {
  get_val <- function(dt, yr, col) {
    v <- dt[yil==yr, get(col)]
    if (!length(v) || is.na(v[1])) return(NA)
    v[1]
  }
  ara_don   <- get_val(don,     y-1, "Ara"); if (is.na(ara_don))   ara_don   <- 0
  oca_don   <- get_val(don,     y,   "Oca"); if (is.na(oca_don))   oca_don   <- 0
  sub_don   <- get_val(don,     y,   "Sub"); if (is.na(sub_don))   sub_don   <- 0
  ara_soguk <- get_val(soguk,   y-1, "Ara"); if (is.na(ara_soguk)) ara_soguk <- 0
  oca_soguk <- get_val(soguk,   y,   "Oca"); if (is.na(oca_soguk)) oca_soguk <- 0
  sub_soguk <- get_val(soguk,   y,   "Sub"); if (is.na(sub_soguk)) sub_soguk <- 0
  ara_min   <- get_val(mintemp, y-1, "Ara")
  oca_min   <- get_val(mintemp, y,   "Oca")
  sub_min   <- get_val(mintemp, y,   "Sub")
  data.table(yil=y,
             kis_don_gun   = ara_don + oca_don + sub_don,
             kis_soguk_gun = ara_soguk + oca_soguk + sub_soguk,
             kis_min_ort   = mean(c(ara_min, oca_min, sub_min), na.rm=TRUE))
}))

# ============================================================
# 2. CLIMATE DATA
# ============================================================
iklim <- as.data.table(read_excel(
  file.path(BASE, "iklim_fst_tam_tablo.xlsx"), sheet=1, skip=2))

setnames(iklim, 1:23, c(
  "sample","ay","yil","onceki_ay",
  "prev_ort","prev_max","prev_min","prev_nem","prev_sicak_gun","prev_soguk_gun",
  "ort_sic","max_sic","min_sic","nem","sicak_gun","soguk_gun","don_gun",
  "delta_ort","delta_max","delta_min","delta_nem","fst","fst_sira"))

iklim <- iklim[!is.na(sample) & sample != "2015_09"]

num_cols <- c("prev_min","prev_max","prev_ort","delta_ort","delta_max","delta_min",
              "min_sic","max_sic","nem","prev_nem","delta_nem","fst")
for (col in num_cols) iklim[[col]] <- suppressWarnings(as.numeric(iklim[[col]]))

iklim[, abs_delta_ort := abs(delta_ort)]
iklim[, abs_delta_max := abs(delta_max)]
iklim[, abs_delta_min := abs(delta_min)]
iklim[, abs_delta_nem := abs(delta_nem)]
iklim[, yil := as.integer(substr(sample, 1, 4))]
iklim <- merge(iklim, kis_skoru, by="yil", all.x=TRUE)

# ============================================================
# 3. DIVERSITY DATA
# ============================================================
samples_list <- c(
  "2015_06","2015_08","2015_10",
  "2016_06","2016_08","2016_09","2016_10",
  "2017_06","2017_10",
  "2018_06","2018_08","2018_09","2018_10",
  "2019_06","2019_08","2019_09","2019_10",
  "2020_06","2020_08","2020_09","2020_10",
  "2021_06","2021_08","2021_09","2021_10"
)

read_div <- function(measure) {
  rbindlist(lapply(samples_list, function(s) {
    f <- file.path(DIR, paste0(s, "_", measure, ".txt"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f, header=FALSE, col.names=c("chr","pos","count","frac","mean"))
    dt <- dt[mean != "na" & !is.na(suppressWarnings(as.numeric(mean)))]
    dt[, mean := as.numeric(mean)]
    dt <- dt[frac >= 0.5]
    if (nrow(dt) == 0) return(NULL)
    data.table(sample=s, mean_val=mean(dt$mean, na.rm=TRUE))
  }), fill=TRUE)
}

div <- merge(read_div("pi")[,.(sample, pi=mean_val)],
             read_div("theta")[,.(sample, theta=mean_val)], by="sample")
div <- merge(div, read_div("D")[,.(sample, D=mean_val)], by="sample")
dat <- merge(div, iklim, by="sample")
cat("Merged rows:", nrow(dat), "\n")

# ============================================================
# 4. CORRELATION PLOT FUNCTION
# ============================================================
make_cor_plot <- function(dat, x_col, y_col, x_lab, y_lab) {
  sub_dat <- dat[!is.na(get(x_col)) & !is.na(get(y_col))]
  if (nrow(sub_dat) < 5) return(NULL)
  ct <- cor.test(sub_dat[[x_col]], sub_dat[[y_col]])
  r  <- ct$estimate
  p  <- ct$p.value
  p_label <- ifelse(p<0.001,"p<0.001",
             ifelse(p<0.01, "p<0.01*",
             ifelse(p<0.05, "p<0.05*", sprintf("p=%.2f",p))))
  col <- ifelse(p<0.05,"red3","grey40")
  ggplot(sub_dat, aes(x=.data[[x_col]], y=.data[[y_col]])) +
    geom_point(aes(color=factor(substr(sample,6,7))), size=3) +
    geom_smooth(method="lm", se=TRUE, color="steelblue", linewidth=0.8) +
    scale_color_manual(values=c("06"="#2ecc71","08"="#e74c3c",
                                "09"="#f39c12","10"="#8e44ad"),
                       labels=c("06"="June","08"="Aug","09"="Sep","10"="Oct"),
                       name="Month") +
    annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5,
             label=sprintf("r=%.2f, %s",r,p_label),
             size=4, color=col, fontface="bold") +
    geom_text(aes(label=sample), size=2.5, vjust=-0.8, color="grey50") +
    labs(title=paste0(y_lab, " vs ", x_lab), x=x_lab, y=y_lab) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(size=9))
}

# ============================================================
# 5. VARIABLE SETS
# ============================================================
prev_vars <- list(
  list(col="prev_min",  lab="Previous Month Min Temp (C)"),
  list(col="prev_max",  lab="Previous Month Max Temp (C)"),
  list(col="delta_max", lab="Delta Max Temp (C)"),
  list(col="delta_min", lab="Delta Min Temp (C)")
)
abs_vars <- list(
  list(col="abs_delta_ort", lab="|Delta| Mean Temp (C)"),
  list(col="abs_delta_max", lab="|Delta| Max Temp (C)"),
  list(col="abs_delta_min", lab="|Delta| Min Temp (C)")
)
winter_vars <- list(
  list(col="kis_don_gun",   lab="Winter Frost Days"),
  list(col="kis_soguk_gun", lab="Winter Cold Days <5C"),
  list(col="kis_min_ort",   lab="Winter Mean Min Temp (C)")
)
current_vars <- list(
  list(col="min_sic", lab="Sampling Month Min Temp (C)"),
  list(col="max_sic", lab="Sampling Month Max Temp (C)"),
  list(col="nem",     lab="Sampling Month Humidity (%)")
)
nem_vars <- list(
  list(col="nem",           lab="Sampling Month Humidity (%)"),
  list(col="prev_nem",      lab="Previous Month Humidity (%)"),
  list(col="delta_nem",     lab="Delta Humidity (%)"),
  list(col="abs_delta_nem", lab="|Delta| Humidity (%)"),
  list(col="delta_ort",     lab="Delta Mean Temp (C)"),
  list(col="abs_delta_ort", lab="|Delta| Mean Temp (C)")
)
div_vars <- list(
  list(col="pi",    lab="Mean nucleotide diversity (pi)"),
  list(col="theta", lab="Mean Watterson's theta"),
  list(col="D",     lab="Mean Tajima's D")
)

make_pdf <- function(filename, climate_vars, ncol_n, width_n=16, height_n=12) {
  plots <- list()
  for (dv in div_vars) for (cv in climate_vars) {
    p <- make_cor_plot(dat, cv$col, dv$col, cv$lab, dv$lab)
    if (!is.null(p)) plots <- c(plots, list(p))
  }
  pdf(file.path(DIR, filename), width=width_n, height=height_n)
  grid.arrange(grobs=plots, ncol=ncol_n)
  dev.off()
  cat("Saved:", filename, "\n")
}

make_pdf("climate_prev_month_vs_diversity.pdf",   prev_vars,    ncol_n=4)
make_pdf("abs_delta_vs_diversity.pdf",             abs_vars,     ncol_n=3, width_n=14)
make_pdf("winter_severity_vs_diversity.pdf",       winter_vars,  ncol_n=3, width_n=14)
make_pdf("current_month_climate_vs_diversity.pdf", current_vars, ncol_n=3, width_n=14)
make_pdf("humidity_delta_vs_diversity.pdf",        nem_vars,     ncol_n=6, width_n=20)

# ============================================================
# 6. CORRELATION TABLE
# ============================================================
all_vars <- c(prev_vars, abs_vars, winter_vars, current_vars, nem_vars)
cor_table <- rbindlist(lapply(div_vars, function(dv) {
  rbindlist(lapply(all_vars, function(cv) {
    sub_dat <- dat[!is.na(get(cv$col)) & !is.na(get(dv$col))]
    if (nrow(sub_dat) < 5) return(NULL)
    ct <- cor.test(sub_dat[[cv$col]], sub_dat[[dv$col]])
    data.table(diversity=dv$col, climate=cv$col,
               r=round(ct$estimate,3), p=round(ct$p.value,4),
               n=nrow(sub_dat), sig=ifelse(ct$p.value<0.05,"*",""))
  }))
}))

cat("\n=== Correlation Table ===\n")
print(cor_table)
fwrite(cor_table, file.path(DIR, "all_climate_diversity_correlations.tsv"), sep="\t")
cat("\n>>> All done!\n")
