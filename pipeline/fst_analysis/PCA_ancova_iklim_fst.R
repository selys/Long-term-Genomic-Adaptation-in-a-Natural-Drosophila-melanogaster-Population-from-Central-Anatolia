library(data.table)
library(ggplot2)
library(ggrepel)
library(dplyr)

# ============================================================
# 1. CLIMATE VERİSİ
# ============================================================
climate <- data.frame(
  Sample = c("2015_06","2015_08","2015_09","2015_10",
            "2016_06","2016_08","2016_09","2016_10",
            "2017_06","2017_10",
            "2018_06","2018_08","2018_09","2018_10",
            "2019_06","2019_08","2019_09","2019_10",
            "2020_06","2020_08","2020_09","2020_10",
            "2021_06","2021_08","2021_09","2021_10"),
  Month = c("June","August","September","October",
         "June","August","September","October",
         "June","October",
         "June","August","September","October",
         "June","August","September","October",
         "June","August","September","October",
         "June","August","September","October"),
  Year = c(2015,2015,2015,2015,
          2016,2016,2016,2016,
          2017,2017,
          2018,2018,2018,2018,
          2019,2019,2019,2019,
          2020,2020,2020,2020,
          2021,2021,2021,2021),
  MeanTemp  = c(18.2,24.6,23.0,14.1, 21.7,25.7,19.0,13.1, 19.7,11.8,
              21.3,25.4,20.3,14.5, 21.5,24.3,20.2,15.5, 20.2,25.2,23.1,16.6,
              19.1,25.3,18.3,12.5),
  MaxSic  = c(25.4,33.2,31.9,21.2, 30.1,34.3,27.7,21.5, 27.8,19.9,
              30.2,34.2,28.8,22.4, 29.9,33.3,29.1,24.7, 28.5,34.4,32.3,25.3,
              27.0,34.2,26.0,20.7),
  MinSic  = c(12.6,16.7,14.9, 8.6, 13.0,17.2,11.4, 6.4, 12.7, 5.7,
              13.6,17.2,12.3, 8.2, 14.7,15.9,12.1, 8.0, 12.8,16.0,14.8, 9.9,
              12.0,16.3,11.4, 5.7),
  Nem     = c(73.8,47.3,43.4,68.4, 53.4,43.3,50.4,56.7, 63.2,61.7,
              59.6,43.9,48.5,68.7, 69.9,49.9,48.9,60.4, 68.5,42.8,51.1,64.0,
              71.0,47.7,61.3,67.3),
  SicakGun = c(0,28,22,0, 15,27,12,0, 8,0, 17,31,10,0, 15,27,15,2, 12,31,22,3,
               10,30, 3,0),
  SogukGun = c(0,0,0,4, 0,0,4,11, 0,17, 0,0,0,6, 0,0,2,9, 0,0,0,0, 3,0,1,14),
  DonGunu  = c(0,0,0,1, 0,0,0,1, 0,1, 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,1)
)

# Season etiketi
climate$Season <- ifelse(climate$Month %in% c("June"), "Summer-Start",
                ifelse(climate$Month %in% c("August","September"), "Summer-Summit",
                "Autumn"))

# ============================================================
# 2. FST VERİSİ — sliding window dosyasından genome-wide ortalama
# ============================================================
cat(">>> FST dosyasi okunuyor...\n")
fst_data <- fread("mega_all_years_26_sliding_fst.txt", data.table = FALSE, header = FALSE)

pop_names <- c("2015_06","2015_08","2015_09","2015_10",
               "2016_06","2016_08","2016_09","2016_10",
               "2017_06","2017_10",
               "2018_06","2018_08","2018_09","2018_10",
               "2019_06","2019_08","2019_09","2019_10",
               "2020_06","2020_08","2020_09","2020_10",
               "2021_06","2021_08","2021_09","2021_10")

sample_cols_idx <- 6:ncol(fst_data)
pair_names <- sub("=.*", "", as.character(fst_data[1, sample_cols_idx]))

fst_numeric <- apply(fst_data[, sample_cols_idx], 2, function(col) {
  as.numeric(sub(".*=", "", as.character(col)))
})

mean_fst <- colMeans(fst_numeric, na.rm = TRUE)
names(mean_fst) <- pair_names

# Her örnek için: diğer tüm örneklerle ortalama FST (mean pairwise FST)
fst_long <- data.frame(
  i = as.integer(sub(":.*", "", names(mean_fst))),
  j = as.integer(sub(".*:", "", names(mean_fst))),
  FST = as.numeric(mean_fst)
)

# Her örnek için ortalama pairwise FST
# fst_long'u iki yönlü mirror yaparak her örneği hem i hem j olarak say
fst_mirror2 <- rbind(
  data.frame(Sample = fst_long$i, FST = fst_long$FST),
  data.frame(Sample = fst_long$j, FST = fst_long$FST)
)
mean_pw_fst <- tapply(fst_mirror2$FST, fst_mirror2$Sample, mean, na.rm = TRUE)
climate$MeanFST <- as.numeric(mean_pw_fst[as.character(1:26)])

cat("FST range:", range(climate$MeanFST, na.rm=TRUE), "\n")

# ============================================================
# 3. PCA
# ============================================================
pca_vars <- c("MeanTemp","MaxSic","MinSic","Nem","SicakGun","SogukGun","DonGunu")
pca_input <- climate[, pca_vars]
pca_res <- prcomp(pca_input, scale. = TRUE)

pca_df <- as.data.frame(pca_res$x[, 1:3])
pca_df$Sample  <- climate$Sample
pca_df$Year    <- factor(climate$Year)
pca_df$Season <- climate$Season
pca_df$MeanFST <- climate$MeanFST

var_exp <- round(summary(pca_res)$importance[2, 1:3] * 100, 1)

# PCA plot — yıla göre renk, Seasone göre şekil
p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Year, shape = Season, label = Sample)) +
  geom_point(aes(size = MeanFST), alpha = 0.85) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  scale_size_continuous(range = c(2, 8), name = "Mean pairwise FST") +
  theme_minimal(base_size = 13) +
  labs(
    title = "PCA: CLIMATE Profile (Yesiloz D. melanogaster)",
    subtitle = "Dot Size = mean pairwise FST",
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme(plot.title = element_text(face = "bold"))

ggsave("PCA_Climate_FST.pdf", p_pca, width = 10, height = 8)
cat(">>> PCA kaydedildi.\n")

# PC yüklemeleri (loadings)
loadings_df <- as.data.frame(pca_res$rotation[, 1:3])
loadings_df$Degisken <- rownames(loadings_df)
cat("\n--- PC Yuklemeleri ---\n")
print(loadings_df)

# ============================================================
# 4. ANCOVA 1: FST ~ Year + CLIMATE değişkenleri
# ============================================================
cat("\n\n============ ANCOVA 1: FST ~ YEAR + CLIMATE ============\n")
ancova_Year <- lm(MeanFST ~ factor(Year) + MeanTemp + Nem + SicakGun + SogukGun,
                 data = climate)
print(summary(ancova_Year))
cat("\n--- ANOVA tablosu ---\n")
print(anova(ancova_Year))

# ============================================================
# 5. ANCOVA 2: FST ~ Season + CLIMATE değişkenleri
# ============================================================
cat("\n\n============ ANCOVA 2: FST ~ SEASON + CLIMATE ============\n")
ancova_Season <- lm(MeanFST ~ Season + MeanTemp + Nem + SicakGun + SogukGun,
                    data = climate)
print(summary(ancova_Season))
cat("\n--- ANOVA tablosu ---\n")
print(anova(ancova_Season))

# ============================================================
# 6. FST ~ PC1 + PC2 regresyon (özet görsel)
# ============================================================
climate$PC1 <- pca_df$PC1
climate$PC2 <- pca_df$PC2

lm_pc <- lm(MeanFST ~ PC1 + PC2, data = climate)
cat("\n\n============ FST ~ PC1 + PC2 ============\n")
print(summary(lm_pc))

p_pc1 <- ggplot(climate, aes(x = PC1, y = MeanFST, color = factor(Year), label = Sample)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "gray40") +
  geom_text_repel(size = 2.8) +
  theme_minimal(base_size = 12) +
  labs(title = "FST ~ PC1 (CLIMATE)", x = paste0("PC1 (", var_exp[1], "%)"),
       y = "Mean pairwise FST", color = "Year")

ggsave("FST_vs_PC1.pdf", p_pc1, width = 9, height = 6)

# ============================================================
# 7. ANCOVA sonuçlarını dosyaya yaz
# ============================================================
sink("ANCOVA_Sonuclari.txt")
cat("========== ANCOVA 1: FST ~ YEAR + CLIMATE ==========\n")
print(summary(ancova_Year))
cat("\nANOVA:\n"); print(anova(ancova_Year))
cat("\n\n========== ANCOVA 2: FST ~ SEASON + CLIMATE ==========\n")
print(summary(ancova_Season))
cat("\nANOVA:\n"); print(anova(ancova_Season))
cat("\n\n========== FST ~ PC1 + PC2 ==========\n")
print(summary(lm_pc))
sink()

cat("\n>>> Tum analizler tamamlandi!\n")
cat("    Ciktilar: PCA_Climate_FST.pdf, FST_vs_PC1.pdf, ANCOVA_Sonuclari.txt\n")


