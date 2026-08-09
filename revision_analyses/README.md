# Revision Analyses

Additional analyses performed in response to reviewer comments on the
BMC Genomics submission (peer review round, 2026).

| Folder | Analysis | Addresses |
|---|---|---|
| `sampling_map/` | Sampling location map and temporal sampling scheme | Reviewer 1: provide a sampling map and a figure illustrating sampling dates |
| `assumption_diagnostics/` | Residual normality diagnostics for the linear models used on allele frequency and F<sub>ST</sub> (Shapiro-Wilk tests), compared against non-parametric alternatives (Spearman, Kruskal-Wallis) | Reviewer 2: it is not clear if the assumptions of the parametric tests are fulfilled |
| `gene_length_control/` | Tests whether gene length confounds raw SNP counts and candidate gene rankings | Reviewer 2: longer genes might harbor more SNPs |
| `go_enrichment_length_adjusted/` | Re-tests GO enrichment results using (a) a length-matched background and (b) a length-adjusted logistic regression across the full gene universe, to determine which enriched terms are robust to gene-length bias | Reviewer 2 (gene length) and Reviewer 1 (caution around GO enrichment as a primary result; Pavlidis et al. 2012) |
| `fst_time_dependence/` | Mantel test (permutation-based, Pearson and Spearman) for whether pairwise F<sub>ST</sub> increases with temporal distance between cohorts | Reviewer 2: the authors state F<sub>ST</sub> changes in a time-dependent manner, but there is no statistical test |

Each folder contains the analysis script(s) (`.sh`/`.R`/`.py`) and the
resulting summary output (`.txt`/`.tsv`). See the manuscript's Methods
and Results sections, and the accompanying point-by-point response, for
how each result is reported and interpreted.

For the core pipeline that produced the manuscript's primary results
(alignment, variant calling, F<sub>ST</sub>, nucleotide diversity, inversion
frequencies, climate ANCOVA, GO enrichment), see `../pipeline/`.
