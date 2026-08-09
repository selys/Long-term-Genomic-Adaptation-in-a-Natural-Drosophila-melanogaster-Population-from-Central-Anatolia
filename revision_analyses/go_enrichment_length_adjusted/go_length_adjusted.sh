#!/bin/bash
#SBATCH -p hamsi
#SBATCH -A ssenkal
#SBATCH -J go_length_adj
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --qos=normal
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/go_length_adj_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/go_length_adj_%j.err

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

LEN_TABLE="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL/HIGH_MODERATE_ONLY/candidate_vs_background_length.tsv"
OUT_DIR="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/GO_ENRICHMENT_length_adjusted"
mkdir -p $OUT_DIR

if [ ! -f "$LEN_TABLE" ]; then
    echo "HATA: $LEN_TABLE yok. Once gene_length_control_HM.sh'i calistir."
    exit 1
fi
echo "Using: $LEN_TABLE (full gene universe, no subsampling)"

cat > $OUT_DIR/go_length_adjusted.py << 'PYEOF'
"""
Definitive length-adjusted GO enrichment test.

Instead of subsampling the background to match candidate-gene lengths
(which controls for length bias but loses statistical power), this uses
the FULL gene universe (all ~14,784 autosomal genes, no genes excluded)
and asks, per GO term, via logistic regression:

    is_candidate ~ in_term + log(gene_length)

If the coefficient on `in_term` is significant AFTER controlling for
log(gene_length), the term is genuinely associated with candidate-gene
status independent of length -- not a length artifact. This keeps full
statistical power while still removing the length confound, unlike the
nearest-neighbor length-matched subsample approach used previously.
"""
import gseapy as gp
import pandas as pd
import numpy as np
from scipy.stats import norm

LEN_TABLE = "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL/HIGH_MODERATE_ONLY/candidate_vs_background_length.tsv"
OUT_DIR = "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/GO_ENRICHMENT_length_adjusted"

df = pd.read_csv(LEN_TABLE, sep="\t")
df = df.dropna(subset=["gene_name", "length_bp"])
df = df[df["gene_name"].apply(lambda x: isinstance(x, str) and len(x) > 0)]
df["log_length"] = np.log(df["length_bp"].astype(float))
df["is_candidate"] = df["is_candidate"].astype(bool).astype(int)

gene_names = df["gene_name"].values
log_length = df["log_length"].values
is_candidate = df["is_candidate"].values
n_genes = len(df)
n_candidates = int(is_candidate.sum())
print(f"Full gene universe: {n_genes} genes ({n_candidates} candidates, no subsampling)")

gene_to_idx = {g: i for i, g in enumerate(gene_names)}


def fit_logit_wald(X, y, ridge=1e-6, max_iter=100, tol=1e-8):
    """Minimal IRLS logistic regression with a Wald test on each coefficient.
    Avoids a statsmodels dependency -- only needs numpy/scipy, both already
    required by this environment for the rest of the pipeline."""
    n, p = X.shape
    beta = np.zeros(p)
    for _ in range(max_iter):
        eta = X @ beta
        eta = np.clip(eta, -30, 30)
        mu = 1.0 / (1.0 + np.exp(-eta))
        W = np.clip(mu * (1 - mu), 1e-8, None)
        XtW = X.T * W
        XtWX = XtW @ X + ridge * np.eye(p)
        z = eta + (y - mu) / W
        XtWz = XtW @ z
        try:
            beta_new = np.linalg.solve(XtWX, XtWz)
        except np.linalg.LinAlgError:
            return None, None
        if np.max(np.abs(beta_new - beta)) < tol:
            beta = beta_new
            break
        beta = beta_new
    eta = np.clip(X @ beta, -30, 30)
    mu = 1.0 / (1.0 + np.exp(-eta))
    W = np.clip(mu * (1 - mu), 1e-8, None)
    XtWX = (X.T * W) @ X + ridge * np.eye(p)
    try:
        cov = np.linalg.inv(XtWX)
    except np.linalg.LinAlgError:
        return beta, None
    se = np.sqrt(np.clip(np.diag(cov), 0, None))
    return beta, se


def bh_correct(pvals):
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.clip(ranked, 0, 1)
    out = np.empty(n)
    out[order] = ranked
    return out


gene_sets = {
    "BP": "GO_Biological_Process_2018",
    "MF": "GO_Molecular_Function_2018",
    "CC": "GO_Cellular_Component_2018",
    "KEGG": "KEGG_2019",
}

all_results = []

for cat, gs_name in gene_sets.items():
    print(f"\n>>> Fetching library: {gs_name} ({cat})")
    library = None
    last_err = None
    for org_str in ["Fly", "fly", "Drosophila"]:
        try:
            library = gp.get_library(name=gs_name, organism=org_str)
            print(f"  Fetched with organism='{org_str}'")
            break
        except Exception as e:
            last_err = e
            continue
    if library is None:
        print(f"  Could not fetch library {gs_name}: {last_err}")
        continue
    print(f"  Terms in library: {len(library)}")

    rows = []
    for term, term_genes in library.items():
        idx = [gene_to_idx[g] for g in term_genes if g in gene_to_idx]
        if len(idx) < 3:
            continue
        in_term = np.zeros(n_genes)
        in_term[idx] = 1
        candidate_hits = int(is_candidate[idx].sum())
        if candidate_hits == 0:
            continue  # nothing to test -- no candidate gene in this term

        X = np.column_stack([np.ones(n_genes), in_term, log_length])
        beta, se = fit_logit_wald(X, is_candidate.astype(float))
        if beta is None or se is None or se[1] == 0 or not np.isfinite(se[1]):
            continue

        z = beta[1] / se[1]
        p = 2 * (1 - norm.cdf(abs(z)))
        rows.append({
            "Gene_set": gs_name, "Ontology": cat, "Term": term,
            "n_term_genes_in_universe": len(idx),
            "candidate_hits": candidate_hits,
            "coef_in_term": beta[1], "se_in_term": se[1],
            "z": z, "P-value": p,
        })

    if not rows:
        print(f"  No testable terms for {cat}")
        continue

    res = pd.DataFrame(rows)
    res["Adjusted P-value"] = bh_correct(res["P-value"].values)
    sig = res[res["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value")
    print(f"  Terms tested: {len(res)} | Significant (q<0.05) after length adjustment: {len(sig)}")

    res.sort_values("P-value").to_csv(f"{OUT_DIR}/GO_{cat}_length_adjusted_ALL_tested.tsv", sep="\t", index=False)
    if len(sig) > 0:
        sig.to_csv(f"{OUT_DIR}/GO_{cat}_length_adjusted_significant.tsv", sep="\t", index=False)
        all_results.append(sig)

if all_results:
    combined = pd.concat(all_results, ignore_index=True).sort_values("Adjusted P-value")
    combined.to_csv(f"{OUT_DIR}/GO_all_significant_length_adjusted.tsv", sep="\t", index=False)
    print(f"\nTotal significant GO terms (length-adjusted, full power): {len(combined)}")
else:
    print("\nNo significant terms found after length adjustment.")

print(f"\n>>> Done: results in {OUT_DIR}")
PYEOF

python $OUT_DIR/go_length_adjusted.py

echo ""
echo ">>> Final verdict table for the 6 terms of interest ==="
cat > $OUT_DIR/final_verdict.py << 'PYEOF'
import pandas as pd
import glob

OUT_DIR = "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/GO_ENRICHMENT_length_adjusted"
TERMS_OF_INTEREST = [
    "integral component of plasma membrane (GO:0005887)",
    "zinc ion binding (GO:0008270)",
    "transition metal ion binding (GO:0046914)",
    "metalloexopeptidase activity (GO:0008235)",
    "metalloaminopeptidase activity (GO:0070006)",
    "carbon-nitrogen ligase activity, with glutamine as amido-N-donor (GO:0016884)",
]

frames = []
for f in glob.glob(f"{OUT_DIR}/GO_*_length_adjusted_ALL_tested.tsv"):
    frames.append(pd.read_csv(f, sep="\t"))
if not frames:
    print("No length-adjusted result files found.")
else:
    all_tested = pd.concat(frames, ignore_index=True)
    hit = all_tested[all_tested["Term"].isin(TERMS_OF_INTEREST)]
    if len(hit) == 0:
        print("None of the terms of interest were testable in the length-adjusted model.")
    else:
        cols = ["Term", "candidate_hits", "coef_in_term", "P-value", "Adjusted P-value"]
        print(hit[cols].sort_values("Adjusted P-value").to_string(index=False))
        print("\nVerdict: Adjusted P-value < 0.05 => survives full-power length adjustment (likely real).")
        print("         Adjusted P-value >= 0.05 => does not survive => treat as a likely length artifact.")
PYEOF

python $OUT_DIR/final_verdict.py

echo ">>> Job finished: $(date)"
