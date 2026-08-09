# Revision Analyses (BMC Genomics resubmission, Aug 2026)

Reviewer 1 & 2 yorumlarına cevaben yapılan ek analizler.

- `gene_length_control/` — Gen uzunluğu ile SNP significance arasındaki
  ilişkinin kontrolü (Reviewer 2: "longer genes might harbor more SNPs").
  Top candidate gene sıralamasının ve GO enrichment sonuçlarının gene
  length bias'ından bağımsız olup olmadığını test eder.
- `go_enrichment_length_adjusted/` — GO enrichment'ın length-matched
  background ve length-adjusted lojistik regresyon ile tekrar test
  edilmesi (Pavlidis et al. 2012 önerisi + Reviewer 2 gene length
  eleştirisi).
- `fst_time_dependence/` — Pairwise FST'nin zamansal mesafeyle
  ilişkisinin Mantel testiyle (permütasyon bazlı, Pearson + Spearman)
  formal olarak test edilmesi (Reviewer 2: "there is no statistical
  test").
- `assumption_diagnostics/` — Parametrik test varsayımlarının
  (normallik, vb.) kontrolü ve parametrik olmayan alternatiflerle
  (Spearman, Kruskal-Wallis) karşılaştırılması (Reviewer 2: "not clear
  if the assumptions of the parametric tests are fulfilled").
- `sampling_map/` — Örnekleme haritası ve zaman çizelgesi figürü
  (Reviewer 1: "provide a sampling map and a small figure illustrating
  the sampling dates").
