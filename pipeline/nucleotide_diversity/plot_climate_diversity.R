library(data.table)
library(ggplot2)
library(gridExtra)
library(readxl)

BASE <- "/arf/scratch/ssenkal/selin/sekanslar"

# ============================================================
# 1. KIŞ SKORU — CSV'lerden oku
# ============================================================
don    <- fread(file.path(BASE, "don_gunler.csv"))
soguk  <- fread(file.path(BASE, "soguk_gunler.csv"))
mintemp <- fread(file.path(BASE, "min_sicaklik.csv"))

# Her yıl için kış skoru: Aralık(t-1) + Ocak(t) + Şubat(t)
kis_skoru <- rbindlist(lapply(2015:2021, function(y) {
  ara_don   <- don[yil == (y-1), Ara];    if (!length(ara_don))   ara_don   <- 0
  oca_don   <- don[yil == y,     Oca];    if (!length(oca_don))   oca_don   <- 0
  sub_don   <- don[yil == y,     Sub];    if (!length(sub_don))   sub_don   <- 0

  ara_soguk <- soguk[yil == (y-1), Ara];  if (!length(ara_soguk)) ara_soguk <- 0
  oca_soguk <- soguk[yil == y,     Oca];  if (!length(oca_soguk)) oca_soguk <- 0
  sub_soguk <- soguk[yil == y,     Sub];  if (!length(sub_soguk)) sub_soguk <- 0

  ara_min   <- mintemp[yil == (y-1), Ara]; if (!length(ara_min))  ara_min   <- NA
  oca_min   <- mintemp[yil == y,     Oca]; if (!length(oca_min))  oca_min   <- NA
  sub_min   <- mintemp[yil == y,     Sub]; if (!length(sub_min))  sub_min   <- NA

  data.table(
    yil           = y,
    kis_don_gun   = ara_don + oca_don + sub_don,
    kis_soguk_gun = ara_soguk + oca_soguk + sub_soguk,
    kis_min_ort   = mean(c(ara_min, oca_min, sub_min), na.rm=TRUE)
  )
}))

cat("Kış skorları:\n"); print(kis_skoru)

# ============================================================
# 2. İKLİM VERİSİ (önceki ay + delta)
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
              "min_sic","max_sic","ort_sic","fst")
for (col in num_cols) iklim[[col]] <- suppressWarnings(as.numeric(iklim[[col]]))
iklim[, yil := as.integer(substr(sample, 1, 4))]
iklim <- merge(iklim, kis_skoru, by="yil", all.x=TRUE)

# ============================================================
# 3. DİVERSİTY VERİSİ
# ============================================================
DIR_1k <- file.path(BASE, "all_years_26/PI_ANALYSIS_1k_ps80")

samples_list <- c("2015_06","2015_08","2015_10",
                  "2016_06","2016_08","2016_09","2016_10",
                  "2017_06","2017_10",
                  "2018_06","2018_08","2018_09","2018_10",
                  "2019_06","2019_08","2019_09","2019_10",
                  "2020_06","2020_08","2020_09","2020_10",
                  "2021_06","2021_08","2021_09","2021_10")

read_div <- function(dir, measure) {
  rbindlist(lapply(samples_list, function(s) {
    f <- file.path(dir, paste0(s, "_", measure, ".txt"))
    if (!file.exists(f)) return(NULL)
    dt <- fread(f, header=FALSE, col.names=c("chr","pos","count","frac","mean"))
    dt <- dt[mean != "na" & !is.na(suppressWarnings(as.numeric(mean)))]
    dt[, mean := as.numeric(mean)]
    dt <- dt[frac >= 0.5]
    if (nrow(dt) == 0) return(NULL)
    data.table(sample=s, mean_val=mean(dt$mean, na.rm=TRUE))
  }), fill=TRUE)
}

div <- merge(read_div(DIR_1k,"pi")[,.(sample,pi=mean_val)],
             read_div(DIR_1k,"theta")[,.(sample,theta=mean_val)], by="sample")
div <- merge(div, read_div(DIR_1k,"D")[,.(sample,D=mean_val)], by="sample")

dat <- merge(div, iklim, by="sample")
cat("\nBirleşik veri:", nrow(dat), "satır\n")

# ============================================================
# 4. KORELASYON PLOT
# ============================================================
make_cor_plot <- function(dat, x_col, y_col, x_lab, y_lab, title) {
  sub_dat <- dat[!is.na(get(x_col)) & !is.na(get(y_col))]
  if (nrow(sub_dat) < 5) return(NULL)
  ct  <- cor.test(sub_dat[[x_col]], sub_dat[[y_col]])
  r   <- ct$estimate
  p   <- ct$p.value
  p_label <- ifelse(p<0.001,"p<0.001",ifelse(p<0.01,"p<0.01",
             ifelse(p<0.05,"p<0.05",sprintf("p=%.2f",p))))
  col <- ifelse(p<0.05,"red3","grey40")
  ggplot(sub_dat, aes(x=.data[[x_col]], y=.data[[y_col]])) +
    geom_point(aes(color=factor(substr(sample,6,7))), size=3) +
    geom_smooth(method="lm", se=TRUE, color="steelblue", linewidth=0.8) +
    scale_color_manual(values=c("6"="#2ecc71","8"="#e74c3c",
                                "9"="#f39c12","10"="#8e44ad"),
                       labels=c("6"="June","8"="Aug","9"="Sep","10"="Oct"),
                       name="Month") +
    annotate("text", x=Inf, y=Inf, hjust=1.1, vjust=1.5,
             label=sprintf("r=%.2f, %s",r,p_label),
             size=4, color=col, fontface="bold") +
    labs(title=title, x=x_lab, y=y_lab) +
    theme_minimal(base_size=10) +
    theme(plot.title=element_text(size=9))
}

# ============================================================
# 5. PLOTLAR
# ============================================================
OUT <- DIR_1k

prev_vars <- list(
  list(col="prev_min",  lab="Prev Month Min Temp (°C)"),
  list(col="prev_max",  lab="Prev Month Max Temp (°C)"),
  list(col="delta_max", lab="Delta Max Temp (°C)"),
  list(col="delta_min", lab="Delta Min Temp (°C)")
)
winter_vars <- list(
  list(col="kis_don_gun",   lab="Winter Frost Days"),
  list(col="kis_soguk_gun", lab="Winter Cold Days <5°C"),
  list(col="kis_min_ort",   lab="Winter Mean Min Temp (°C)")
)
div_vars <- list(
  list(col="pi",    lab="Mean pi"),
  list(col="theta", lab="Mean theta"),
  list(col="D",     lab="Mean Tajima's D")
)

# PDF 1: önceki ay iklimi
pdf(file.path(OUT,"climate_prev_month_vs_diversity.pdf"), width=16, height=12)
plots <- list()
for (dv in div_vars) for (cv in prev_vars) {
  p <- make_cor_plot(dat, cv$col, dv$col, cv$lab, dv$lab,
                     paste0(dv$lab,"\nvs ",cv$lab))
  if (!is.null(p)) plots <- c(plots, list(p))
}
grid.arrange(grobs=plots, ncol=4)
dev.off()
cat("Saved: climate_prev_month_vs_diversity.pdf\n")

# PDF 2: kış sertliği
pdf(file.path(OUT,"winter_severity_vs_diversity.pdf"), width=14, height=12)
plots2 <- list()
for (dv in div_vars) for (cv in winter_vars) {
  p <- make_cor_plot(dat, cv$col, dv$col, cv$lab, dv$lab,
                     paste0(dv$lab,"\nvs ",cv$lab))
  if (!is.null(p)) plots2 <- c(plots2, list(p))
}
grid.arrange(grobs=plots2, ncol=3)
dev.off()
cat("Saved: winter_severity_vs_diversity.pdf\n")

# Korelasyon tablosu
all_vars <- c(prev_vars, winter_vars)
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

cat("\n=== Korelasyon Tablosu ===\n"); print(cor_table)
fwrite(cor_table, file.path(OUT,"climate_diversity_correlations.tsv"), sep="\t")
cat("\nDone!\n")
