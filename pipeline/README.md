# Core Analysis Pipeline

This folder contains the scripts that produced ALL of the manuscript's
primary results (for additional control analyses added during revision,
see `../revision_analyses/`).

- `variant_calling/` — mpileup, VarScan, SnpEff, sync file creation
- `fst_analysis/` — PoPoolation2 sliding-window F<sub>ST</sub>, F<sub>ST</sub> heatmap, climate PCA
- `nucleotide_diversity/` — π, Watterson's θ, and Tajima's D via PoPoolation1
  (pool-size 80/160, 1 kb sliding windows, autosomal chromosomes only)
- `inversion_analysis/` — Frequency estimation for 7 cosmopolitan
  chromosomal inversions using the Kapun et al. (2014) diagnostic
  marker panel
- `climate_ancova/` — SNP-level allele frequency ~ season + year +
  climate PCA models, joint multiple-testing correction, Manhattan plot
- `go_enrichment/` — Original (whole-genome background) GO enrichment

If a script could not be located automatically when this folder was
built (flagged `[MISSING]` in the build log), its real location on the
compute cluster still needs to be found and copied in manually.
