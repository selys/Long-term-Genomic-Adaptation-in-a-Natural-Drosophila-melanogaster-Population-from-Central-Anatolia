# Core Analysis Pipeline

Bu klasör, manuscript'in TÜM ana bulgularını üreten script'leri içerir
(revizyon sırasında eklenen ek kontrol analizleri için bkz. `../revision_analyses/`).

- `variant_calling/` — mpileup, VarScan, SnpEff, sync dosyası oluşturma
- `fst_analysis/` — PoPoolation2 sliding-window FST, FST heatmap, climate PCA
- `nucleotide_diversity/` — PoPoolation1 ile π, Watterson's θ, Tajima's D
  (pool-size 80/160, 1kb pencere, otozomal kromozomlar)
- `inversion_analysis/` — Kapun et al. (2014) diagnostic marker paneli ile
  7 kozmopolit inversion'ın frekans tahmini
- `climate_ancova/` — SNP-bazlı allele frequency ~ season + year + climate
  PCA modelleri, joint multiple-testing correction, Manhattan plot
- `go_enrichment/` — Orijinal (whole-genome background) GO enrichment

Bazı dosyalar otomatik bulunamadıysa (`[YOK]` uyarısı), ARF'taki gerçek
konumlarını bulup elle bu klasörlere kopyalaman gerekir.
