# =============================================================================
# SNP Allele Frequency - Climate ANCOVA Analysis
# Models: Season-only, Year-only, Season*Year, Season+Year+Climate(PCA)
# Input:  mega_all_years_freq_diff_rc
# Output: CLIMATE_ANALYSIS/
# =============================================================================

library(data.table)
library(foreach)
library(doParallel)
library(car)

# -----------------------------------------------------------------------------
# 0. SETUP
# -----------------------------------------------------------------------------
registerDoParallel(cores = 40)

INPUT_FILE <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/MAF_ANALYSIS/mega_all_years_freq_diff_rc"
OUT_DIR    <- "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# -----------------------------------------------------------------------------
# 1. SAMPLE METADATA (26 samples, 2017 has only June + October)
# -----------------------------------------------------------------------------
sample_meta <- data.table(
  sample_id = paste0("maa_", 1:26),
  year      = c(rep(2015,4), rep(2016,4), 2017, 2017,
                rep(2018,4), rep(2019,4), rep(2020,4), rep(2021,4)),
  month     = c(6,8,9,10,  6,8,9,10,  6,10,
                6,8,9,10,  6,8,9,10,  6,8,9,10,  6,8,9,10)
)
sample_meta[, season := factor(month, levels=c(6,8,9,10),
                               labels=c("June","August","September","October"))]
sample_meta[, year_val := as.numeric(year)]

cat("Sample metadata:\n")
print(sample_meta)

# -----------------------------------------------------------------------------
# 2. CLIMATE DATA
# Months in dataset: 6=June, 8=August, 9=September, 10=October
# 2017 has no August/September → those rows get NA for climate vars
# -----------------------------------------------------------------------------

# Helper: build a (year x month) lookup table
make_clim <- function(varname, vals_list) {
  rows <- lapply(names(vals_list), function(yr) {
    m <- vals_list[[yr]]
    data.table(year=as.integer(yr), month=as.integer(names(m)),
               value=unlist(m, use.names=FALSE))
  })
  dt <- rbindlist(rows)
  setnames(dt, "value", varname)
  dt
}

avg_temp <- make_clim("avg_temp", list(
  "2015"=c("6"=18.2,"8"=24.6,"9"=23.0,"10"=14.1),
  "2016"=c("6"=21.7,"8"=25.7,"9"=19.0,"10"=13.1),
  "2017"=c("6"=19.7,         "10"=11.8),
  "2018"=c("6"=21.3,"8"=25.4,"9"=20.3,"10"=14.5),
  "2019"=c("6"=21.5,"8"=24.3,"9"=20.2,"10"=15.5),
  "2020"=c("6"=20.2,"8"=25.2,"9"=23.1,"10"=16.6),
  "2021"=c("6"=19.1,"8"=25.3,"9"=18.3,"10"=12.5)))

min_temp <- make_clim("min_temp", list(
  "2015"=c("6"=12.6,"8"=16.7,"9"=14.9,"10"=8.6),
  "2016"=c("6"=13.0,"8"=17.2,"9"=11.4,"10"=6.4),
  "2017"=c("6"=12.7,         "10"=5.7),
  "2018"=c("6"=13.6,"8"=17.2,"9"=12.3,"10"=8.2),
  "2019"=c("6"=14.7,"8"=15.9,"9"=12.1,"10"=8.0),
  "2020"=c("6"=12.8,"8"=16.0,"9"=14.8,"10"=9.9),
  "2021"=c("6"=12.0,"8"=16.3,"9"=11.4,"10"=5.7)))

max_temp <- make_clim("max_temp", list(
  "2015"=c("6"=25.4,"8"=33.2,"9"=31.9,"10"=21.2),
  "2016"=c("6"=30.1,"8"=34.3,"9"=27.7,"10"=21.5),
  "2017"=c("6"=27.8,         "10"=19.9),
  "2018"=c("6"=30.2,"8"=34.2,"9"=28.8,"10"=22.4),
  "2019"=c("6"=29.9,"8"=33.3,"9"=29.1,"10"=24.7),
  "2020"=c("6"=28.5,"8"=34.4,"9"=32.3,"10"=25.3),
  "2021"=c("6"=27.0,"8"=34.2,"9"=26.0,"10"=20.7)))

precip <- make_clim("precipitation", list(
  "2015"=c("6"=152.6,"8"=21.2,"9"=2.7, "10"=21.7),
  "2016"=c("6"=16.1, "8"=4.8, "9"=32.3,"10"=2.4),
  "2017"=c("6"=91.2,          "10"=33.4),
  "2018"=c("6"=22.3, "8"=52.6,"9"=2.0, "10"=67.2),
  "2019"=c("6"=89.6, "8"=9.3, "9"=1.0, "10"=13.8),
  "2020"=c("6"=62.5, "8"=0.4, "9"=8.9, "10"=39.2),
  "2021"=c("6"=118.8,"8"=4.4, "9"=10.4,"10"=25.5)))

hot_days <- make_clim("hot_days", list(
  "2015"=c("6"=0, "8"=28,"9"=22,"10"=0),
  "2016"=c("6"=15,"8"=27,"9"=12,"10"=0),
  "2017"=c("6"=8,         "10"=0),
  "2018"=c("6"=17,"8"=31,"9"=10,"10"=0),
  "2019"=c("6"=15,"8"=27,"9"=15,"10"=2),
  "2020"=c("6"=12,"8"=31,"9"=22,"10"=3),
  "2021"=c("6"=10,"8"=30,"9"=3, "10"=0)))

avg_humidity <- make_clim("avg_humidity", list(
  "2015"=c("6"=73.8,"8"=47.3,"9"=43.4,"10"=68.4),
  "2016"=c("6"=53.4,"8"=43.3,"9"=50.4,"10"=56.7),
  "2017"=c("6"=63.2,         "10"=61.7),
  "2018"=c("6"=59.6,"8"=43.9,"9"=48.5,"10"=68.7),
  "2019"=c("6"=69.9,"8"=49.9,"9"=48.9,"10"=60.4),
  "2020"=c("6"=68.5,"8"=42.8,"9"=51.1,"10"=64.0),
  "2021"=c("6"=71.0,"8"=47.7,"9"=61.3,"10"=67.3)))

min_humidity <- make_clim("min_humidity", list(
  "2015"=c("6"=38.9,"8"=20.4,"9"=18.5,"10"=38.5),
  "2016"=c("6"=22.8,"8"=17.6,"9"=22.6,"10"=29.0),
  "2017"=c("6"=32.5,         "10"=31.7),
  "2018"=c("6"=27.3,"8"=18.5,"9"=23.9,"10"=37.7),
  "2019"=c("6"=34.5,"8"=23.8,"9"=23.3,"10"=31.0),
  "2020"=c("6"=33.7,"8"=20.3,"9"=24.6,"10"=36.3),
  "2021"=c("6"=38.7,"8"=21.8,"9"=35.9,"10"=41.3)))

max_humidity <- make_clim("max_humidity", list(
  "2015"=c("6"=98.9,"8"=77.2,"9"=71.2,"10"=91.9),
  "2016"=c("6"=89.1,"8"=73.2,"9"=80.6,"10"=82.8),
  "2017"=c("6"=89.0,         "10"=85.3),
  "2018"=c("6"=89.8,"8"=72.1,"9"=75.9,"10"=91.4),
  "2019"=c("6"=96.6,"8"=78.0,"9"=76.4,"10"=88.1),
  "2020"=c("6"=94.9,"8"=69.9,"9"=80.6,"10"=86.4),
  "2021"=c("6"=97.7,"8"=76.0,"9"=87.2,"10"=91.4)))

# Merge all climate tables
clim_all <- Reduce(function(a,b) merge(a,b,by=c("year","month"),all=TRUE),
                   list(avg_temp, min_temp, max_temp, precip,
                        hot_days, avg_humidity, min_humidity, max_humidity))

# Derived variables
clim_all[, temp_amplitude   := max_temp - min_temp]
clim_all[, humidity_range   := max_humidity - min_humidity]

# Year-to-year delta (previous year same month)
setorder(clim_all, month, year)
clim_all[, delta_avg_temp   := avg_temp   - shift(avg_temp),   by=month]
clim_all[, delta_min_temp   := min_temp   - shift(min_temp),   by=month]
clim_all[, delta_max_temp   := max_temp   - shift(max_temp),   by=month]
clim_all[, delta_humidity   := avg_humidity - shift(avg_humidity), by=month]
clim_all[, delta_precip     := precipitation - shift(precipitation), by=month]

CLIM_VARS <- c("avg_temp","min_temp","max_temp","precipitation",
               "hot_days","avg_humidity","min_humidity","max_humidity",
               "temp_amplitude","humidity_range",
               "delta_avg_temp","delta_min_temp","delta_max_temp",
               "delta_humidity","delta_precip")

# Merge climate into sample metadata
sample_meta <- merge(sample_meta, clim_all[, c("year","month",CLIM_VARS), with=FALSE],
                     by=c("year","month"), all.x=TRUE)

# -----------------------------------------------------------------------------
# 3. PCA ON CLIMATE VARIABLES (across the 26 samples)
# -----------------------------------------------------------------------------
clim_matrix <- as.matrix(sample_meta[, CLIM_VARS, with=FALSE])
# Use only complete rows for PCA, then predict for all
complete_rows <- complete.cases(clim_matrix)
pca_fit <- prcomp(clim_matrix[complete_rows,], scale.=TRUE, center=TRUE)

cat("\nPCA variance explained:\n")
pct <- summary(pca_fit)$importance[2,] * 100
print(round(pct[1:5], 1))

# Predict PC scores for all 26 samples (NA rows get NA scores)
pc_scores <- predict(pca_fit, newdata=clim_matrix)
sample_meta[, PC1 := pc_scores[,1]]
sample_meta[, PC2 := pc_scores[,2]]
sample_meta[, PC3 := pc_scores[,3]]

cat("\nPC1+PC2 variance explained:", round(sum(pct[1:2]),1), "%\n")

# -----------------------------------------------------------------------------
# 4. READ INPUT FILE
# -----------------------------------------------------------------------------
cat("\n>>> Reading input file...\n")

# Header line starts with ##; fread skips it as a comment.
# Read column names manually, then pass directly to fread.
raw_header <- readLines(INPUT_FILE, n = 1)
col_names  <- strsplit(sub("^##", "", raw_header), "\t")[[1]]
cat("Columns detected:", length(col_names), "\n")

# skip=1 skips the ## header line; we supply col.names ourselves
dt <- fread(INPUT_FILE, skip = 1, header = FALSE, col.names = col_names)

# Parse allele frequencies for each sample
cat(">>> Parsing allele frequencies...\n")
y_list <- vector("list", nrow(sample_meta))

for (i in seq_len(nrow(sample_meta))) {
  samp  <- sample_meta$sample_id[i]
  if (!samp %in% names(dt)) next
  counts <- tstrsplit(dt[[samp]], "/", type.convert=TRUE)
  tot    <- counts[[2]]
  af     <- round(counts[[1]] / tot, 6)
  
  tmp <- dt[, .(chr, pos)]
  tmp[, AF       := af]
  tmp[, tot_c    := tot]
  tmp[, season   := sample_meta$season[i]]
  tmp[, year_val := sample_meta$year_val[i]]
  tmp[, PC1      := sample_meta$PC1[i]]
  tmp[, PC2      := sample_meta$PC2[i]]
  tmp[, PC3      := sample_meta$PC3[i]]
  tmp[, avg_temp     := sample_meta$avg_temp[i]]
  tmp[, avg_humidity := sample_meta$avg_humidity[i]]
  tmp[, temp_amplitude := sample_meta$temp_amplitude[i]]
  
  y_list[[i]] <- tmp[tot_c >= 15]
  rm(tmp, counts, af, tot); gc()
}

rm(dt); gc()
y_master <- rbindlist(y_list)
rm(y_list); gc()

snps <- unique(y_master[, .(chr, pos)])
cat("Total SNPs:", nrow(snps), "\n")

# -----------------------------------------------------------------------------
# 5. PARALLEL MODEL FITTING
# -----------------------------------------------------------------------------
cat("\n>>> Running parallel ANCOVA models...\n")

res <- foreach(i = 1:nrow(snps),
               .combine  = rbind,
               .packages = c("data.table","car"),
               .errorhandling = "remove") %dopar% {

  s   <- snps[i]
  sub <- y_master[chr == s$chr & pos == s$pos]
  if (nrow(sub) < 10) return(NULL)

  tryCatch({
    # ── Model 1: Season only ──────────────────────────────────────────────
    m1   <- lm(AF ~ season, data=sub)
    p1   <- anova(m1)$`Pr(>F)`[1]

    # ── Model 2: Year only ───────────────────────────────────────────────
    m2   <- lm(AF ~ year_val, data=sub)
    p2   <- anova(m2)$`Pr(>F)`[1]

    # ── Model 3: Season * Year (ANCOVA) ──────────────────────────────────
    m3   <- lm(AF ~ season * year_val, data=sub)
    a3   <- Anova(m3, type="III")
    p3_season <- a3$`Pr(>F)`[2]   # season
    p3_year   <- a3$`Pr(>F)`[3]   # year_val
    p3_inter  <- a3$`Pr(>F)`[4]   # season:year_val

    # ── Model 4: Season + Year + Climate PCA ─────────────────────────────
    # Use PC1+PC2; drop rows where PCs are NA (2017 Aug/Sep)
    sub4 <- sub[!is.na(PC1) & !is.na(PC2)]
    if (nrow(sub4) < 10) {
      p4_season <- NA_real_; p4_year <- NA_real_
      p4_PC1 <- NA_real_;    p4_PC2 <- NA_real_
    } else {
      m4   <- lm(AF ~ season + year_val + PC1 + PC2, data=sub4)
      a4   <- Anova(m4, type="III")
      p4_season <- a4$`Pr(>F)`[2]
      p4_year   <- a4$`Pr(>F)`[3]
      p4_PC1    <- a4$`Pr(>F)`[4]
      p4_PC2    <- a4$`Pr(>F)`[5]
    }

    data.table(
      chr=s$chr, pos=s$pos,
      p1_season=p1,
      p2_year=p2,
      p3_season=p3_season, p3_year=p3_year, p3_interaction=p3_inter,
      p4_season=p4_season, p4_year=p4_year,
      p4_PC1=p4_PC1, p4_PC2=p4_PC2
    )
  }, error=function(e) NULL)
}

cat("Models done. Applying BH correction...\n")

# -----------------------------------------------------------------------------
# 6. BH CORRECTION (separately per test column)
# -----------------------------------------------------------------------------
p_cols <- c("p1_season","p2_year",
            "p3_season","p3_year","p3_interaction",
            "p4_season","p4_year","p4_PC1","p4_PC2")

for (pc in p_cols) {
  qc <- sub("^p", "q", pc)
  res[, (qc) := p.adjust(get(pc), method="BH")]
}

# -----------------------------------------------------------------------------
# 7. SAVE FULL RESULTS
# -----------------------------------------------------------------------------
full_out <- file.path(OUT_DIR, "snp_climate_ancova_full.tsv")
fwrite(res, full_out, sep="\t")
cat("Full results saved:", full_out, "\n")

# -----------------------------------------------------------------------------
# 8. SIGNIFICANT SNP SUBSETS (q < 0.05 per test)
# -----------------------------------------------------------------------------
sig_tests <- list(
  season_only      = res[q1_season      < 0.05, .(chr, pos, p1_season, q1_season)],
  year_only        = res[q2_year        < 0.05, .(chr, pos, p2_year,   q2_year)],
  ancova_season    = res[q3_season      < 0.05, .(chr, pos, p3_season, q3_season,
                                                   p3_year, q3_year,
                                                   p3_interaction, q3_interaction)],
  ancova_year      = res[q3_year        < 0.05, .(chr, pos, p3_season, q3_season,
                                                   p3_year, q3_year,
                                                   p3_interaction, q3_interaction)],
  ancova_interact  = res[q3_interaction < 0.05, .(chr, pos, p3_season, q3_season,
                                                   p3_year, q3_year,
                                                   p3_interaction, q3_interaction)],
  climate_PC1      = res[q4_PC1         < 0.05, .(chr, pos, p4_season, q4_season,
                                                   p4_year, q4_year,
                                                   p4_PC1, q4_PC1,
                                                   p4_PC2, q4_PC2)],
  climate_PC2      = res[q4_PC2         < 0.05, .(chr, pos, p4_season, q4_season,
                                                   p4_year, q4_year,
                                                   p4_PC1, q4_PC1,
                                                   p4_PC2, q4_PC2)]
)

for (nm in names(sig_tests)) {
  out_path <- file.path(OUT_DIR, paste0("sig_", nm, ".tsv"))
  fwrite(sig_tests[[nm]], out_path, sep="\t")
  cat(sprintf("  %-20s → %d significant SNPs → %s\n",
              nm, nrow(sig_tests[[nm]]), out_path))
}

# -----------------------------------------------------------------------------
# 9. PCA LOADINGS REPORT
# -----------------------------------------------------------------------------
loadings_out <- file.path(OUT_DIR, "climate_pca_loadings.tsv")
loadings_dt  <- as.data.table(pca_fit$rotation[, 1:3], keep.rownames=TRUE)
setnames(loadings_dt, "rn", "climate_variable")
fwrite(loadings_dt, loadings_out, sep="\t")
cat("PCA loadings saved:", loadings_out, "\n")

cat("\n>>> All done:", format(Sys.time()), "\n")
