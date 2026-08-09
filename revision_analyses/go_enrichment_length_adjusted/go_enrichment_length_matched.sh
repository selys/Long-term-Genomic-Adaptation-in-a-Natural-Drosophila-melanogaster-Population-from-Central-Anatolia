#!/bin/bash
#SBATCH -p hamsi
#SBATCH -A ssenkal
#SBATCH -J go_length_matched
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --mem=16G
#SBATCH -t 02:00:00
#SBATCH --qos=normal
#SBATCH -o /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/go_length_matched_%j.out
#SBATCH -e /arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/go_length_matched_%j.err

source /arf/scratch/ssenkal/selin/programs/miniconda3/etc/profile.d/conda.sh
conda activate drosophila_env

# -----------------------------------------------------------------------
# NEW output folder -- existing GO_ENRICHMENT_* results are untouched.
# -----------------------------------------------------------------------
GLC_DIR="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL"
LEN_TABLE="$GLC_DIR/HIGH_MODERATE_ONLY/candidate_vs_background_length.tsv"
OUT_DIR="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/GO_ENRICHMENT_length_matched"
mkdir -p $OUT_DIR

if [ ! -f "$LEN_TABLE" ]; then
    echo "HATA: $LEN_TABLE yok. Once gene_length_control_HM.sh'i calistir."
    exit 1
fi
echo "Using length table: $LEN_TABLE"

cat > $OUT_DIR/build_matched_background.py << 'PYEOF'
"""
Build a length-matched background gene set for GO enrichment.
For each of the 395 candidate genes, finds the nearest non-candidate genes
by log(gene length) (greedy nearest-neighbor, without replacement where
possible) to form a control pool with a similar length distribution.
This directly answers: "is the enrichment still there once genes are
compared only against other genes of similar size?"
"""
import pandas as pd
import numpy as np

IN_FILE = "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/GENE_LENGTH_CONTROL/HIGH_MODERATE_ONLY/candidate_vs_background_length.tsv"
OUT_DIR = "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/GO_ENRICHMENT_length_matched"
N_MATCH_PER_CANDIDATE = 5   # background pool size = ~5x candidate set

df = pd.read_csv(IN_FILE, sep="\t")

# Drop rows with missing/non-string gene names BEFORE anything else -- a
# handful of GTF records apparently had unresolved names that came through
# as NaN, which crashed the file-write step downstream.
n_before = len(df)
df = df.dropna(subset=["gene_name", "length_bp"])
df = df[df["gene_name"].apply(lambda x: isinstance(x, str) and len(x) > 0)]
n_after = len(df)
print(f"Dropped {n_before - n_after} rows with missing/invalid gene_name (kept {n_after})")

df["log_length"] = np.log(df["length_bp"])

candidates = df[df["is_candidate"] == True].copy()
pool = df[df["is_candidate"] == False].copy().reset_index(drop=True)

print(f"Candidate genes: {len(candidates)}")
print(f"Background pool available: {len(pool)}")

pool_lengths = pool["log_length"].values
pool_used = np.zeros(len(pool), dtype=bool)
matched_genes = []

rng = np.random.default_rng(42)
order = rng.permutation(len(candidates))  # randomize match order to avoid systematic bias

for idx in order:
    row = candidates.iloc[idx]
    target = row["log_length"]
    # distance to every pool gene, mask out already-used ones heavily
    dist = np.abs(pool_lengths - target)
    dist_masked = np.where(pool_used, np.inf, dist)
    # take N_MATCH_PER_CANDIDATE nearest available genes
    nearest_idx = np.argsort(dist_masked)[:N_MATCH_PER_CANDIDATE]
    for i in nearest_idx:
        if not pool_used[i]:
            pool_used[i] = True
            matched_genes.append(pool.iloc[i]["gene_name"])

matched_genes = list(dict.fromkeys(matched_genes))  # dedupe, preserve order
print(f"Matched background genes selected: {len(matched_genes)}")

candidate_list = candidates["gene_name"].tolist()
background_list = list(dict.fromkeys(candidate_list + matched_genes))  # Enrichr background must include gene_list itself

# Final safety net: force everything to str and drop any stray empties
# before writing, so a single bad value can't crash the write again.
candidate_list = [str(g) for g in candidate_list if pd.notna(g) and str(g).strip()]
background_list = [str(g) for g in background_list if pd.notna(g) and str(g).strip()]

if len(background_list) < len(candidate_list):
    raise RuntimeError(
        f"Length-matched background ({len(background_list)} genes) is smaller than "
        f"the candidate set ({len(candidate_list)} genes) -- something went wrong "
        f"in matching. Aborting before writing any files."
    )

with open(f"{OUT_DIR}/candidate_gene_list.txt", "w") as f:
    f.write("\n".join(candidate_list) + "\n")
with open(f"{OUT_DIR}/length_matched_background.txt", "w") as f:
    f.write("\n".join(background_list) + "\n")

print(f"Wrote {len(candidate_list)} candidate genes and {len(background_list)} background genes to disk.")

# Sanity check: compare length distributions of candidates vs matched background
import scipy.stats as stats
matched_lengths = pool[pool["gene_name"].isin(matched_genes)]["length_bp"]
candidate_lengths = candidates["length_bp"]
w, p = stats.mannwhitneyu(candidate_lengths, matched_lengths, alternative="two-sided")
print(f"\nSanity check -- candidate vs MATCHED background length (should now be non-significant):")
print(f"  Median candidate length: {candidate_lengths.median():.0f} bp")
print(f"  Median matched background length: {matched_lengths.median():.0f} bp")
print(f"  Mann-Whitney U={w:.0f}, p={p:.4g}")
PYEOF

echo ">>> Building length-matched background..."
python $OUT_DIR/build_matched_background.py

BG_LINES=$(wc -l < "$OUT_DIR/length_matched_background.txt" 2>/dev/null || echo 0)
CAND_LINES=$(wc -l < "$OUT_DIR/candidate_gene_list.txt" 2>/dev/null || echo 0)
echo "Background file line count: $BG_LINES / Candidate file line count: $CAND_LINES"
if [ "$BG_LINES" -lt "$CAND_LINES" ]; then
    echo "HATA: length_matched_background.txt beklenenden kucuk/bos ($BG_LINES satir)."
    echo "GO enrichment'a gecmiyorum -- once build_matched_background.py'nin hata vermedigini dogrula."
    exit 1
fi

echo ">>> Cleaning any stale results from previous (broken) runs..."
rm -f "$OUT_DIR"/GO_*_length_matched_results.tsv "$OUT_DIR"/GO_all_significant_length_matched.tsv

echo ">>> Running GO enrichment with length-matched background (manual Fisher's exact test)..."
cat > $OUT_DIR/go_enrichment_length_matched.py << 'PYEOF'
"""
GO enrichment for the 395 HIGH/MODERATE candidate genes using a
length-matched custom background, computed MANUALLY via Fisher's exact
test rather than relying on gseapy's enrichr(background=...) parameter
(which errored out with a KeyError in this environment/version -- likely
an API mismatch for list-based custom backgrounds).

Logic per GO term:
    universe = length-matched background genes (INCLUDES the candidates
               themselves, since Enrichr-style background conventions
               treat the gene list as a subset of the background)
        a = candidate genes annotated with this term
        b = candidate genes NOT annotated with this term
        c = background-only genes annotated with this term
        d = background-only genes NOT annotated with this term
    Fisher's exact test (two-sided) on [[a,b],[c,d]], then BH correction
    across all terms that had at least one candidate-gene overlap.
"""
import gseapy as gp
import pandas as pd
import numpy as np
from scipy.stats import fisher_exact

def bh_correct(pvals):
    """Manual Benjamini-Hochberg FDR correction (avoids a statsmodels dependency)."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order] * n / (np.arange(n) + 1)
    # enforce monotonicity (standard BH step-up)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.clip(ranked, 0, 1)
    out = np.empty(n)
    out[order] = ranked
    return out

OUT_DIR = "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/GO_ENRICHMENT_length_matched"

with open(f"{OUT_DIR}/candidate_gene_list.txt") as f:
    candidates = set(line.strip() for line in f if line.strip())
with open(f"{OUT_DIR}/length_matched_background.txt") as f:
    background = set(line.strip() for line in f if line.strip())

# background must be a superset of candidates (Enrichr convention)
background = background | candidates
background_only = background - candidates

print(f"Candidate genes: {len(candidates)}")
print(f"Full background (incl. candidates): {len(background)}")
print(f"Background-only (non-candidate) genes: {len(background_only)}")

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
    for org_str in ["Fly", "fly", "Drosophila"]:
        try:
            library = gp.get_library(name=gs_name, organism=org_str)
            print(f"  Fetched with organism='{org_str}'")
            break
        except Exception as e:
            last_err = e
            continue
    if library is None:
        print(f"  Could not fetch library {gs_name} with any organism string tried: {last_err}")
        continue
    print(f"  Terms in library: {len(library)}")

    rows = []
    for term, term_genes in library.items():
        term_genes = set(term_genes)
        # restrict the term's genes to our background universe only --
        # genes not in our background shouldn't count either direction
        term_in_bg = term_genes & background
        if len(term_in_bg) == 0:
            continue  # term has no representation in our background at all

        a = len(term_in_bg & candidates)
        if a == 0:
            continue  # no candidate genes hit this term -- not enrichable
        b = len(candidates) - a
        c = len(term_in_bg & background_only)
        d = len(background_only) - c

        odds_ratio, p = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append({
            "Gene_set": gs_name, "Ontology": cat, "Term": term,
            "Overlap": f"{a}/{len(term_in_bg)}",
            "Candidate_hits": a, "Background_only_hits": c,
            "Odds_Ratio": odds_ratio, "P-value": p,
        })

    if not rows:
        print(f"  No terms with any candidate-gene overlap for {cat}")
        continue

    res = pd.DataFrame(rows)
    res["Adjusted P-value"] = bh_correct(res["P-value"].values)
    sig = res[res["Adjusted P-value"] < 0.05].sort_values("Adjusted P-value")
    print(f"  Terms tested: {len(res)} | Significant (q<0.05): {len(sig)}")

    res.sort_values("P-value").to_csv(f"{OUT_DIR}/GO_{cat}_length_matched_ALL_tested.tsv", sep="\t", index=False)
    if len(sig) > 0:
        sig.to_csv(f"{OUT_DIR}/GO_{cat}_length_matched_results.tsv", sep="\t", index=False)
        all_results.append(sig)

if all_results:
    combined = pd.concat(all_results, ignore_index=True).sort_values("Adjusted P-value")
    combined.to_csv(f"{OUT_DIR}/GO_all_significant_length_matched.tsv", sep="\t", index=False)
    print(f"\nTotal significant GO terms (length-matched background): {len(combined)}")
else:
    print("\nNo significant terms found with length-matched background.")

print(f"\n>>> Done: results in {OUT_DIR}")
PYEOF

python $OUT_DIR/go_enrichment_length_matched.py

echo ">>> Comparing against the ORIGINAL (whole-genome background) GO results..."
ORIGINAL_GO="/arf/scratch/ssenkal/selin/sekanslar/all_years_26/GO_ENRICHMENT/GO_all_significant.tsv"
MATCHED_GO="$OUT_DIR/GO_all_significant_length_matched.tsv"

cat > $OUT_DIR/compare_terms.py << 'PYEOF'
import pandas as pd
import sys

orig_path = sys.argv[1]
matched_path = sys.argv[2]

try:
    orig = pd.read_csv(orig_path, sep="\t")
    orig_terms = set(orig["Term"])
except Exception as e:
    print(f"Could not read original results ({e}); skipping comparison.")
    sys.exit(0)

print("Original (whole-genome background) significant terms:")
for t in sorted(orig_terms):
    print(f"  - {t}")

try:
    matched = pd.read_csv(matched_path, sep="\t")
    matched_terms = set(matched["Term"])
except FileNotFoundError:
    matched_terms = set()

print("\nLength-matched-background significant terms:")
if matched_terms:
    for t in sorted(matched_terms):
        print(f"  - {t}")
else:
    print("  (none)")

survived = orig_terms & matched_terms
lost = orig_terms - matched_terms
print(f"\n>>> SURVIVED length-matched background ({len(survived)}):")
for t in sorted(survived):
    print(f"  [OK]   {t}")
print(f"\n>>> LOST after length-matched background ({len(lost)}) -- likely length artifacts:")
for t in sorted(lost):
    print(f"  [DROP] {t}")
PYEOF

python $OUT_DIR/compare_terms.py "$ORIGINAL_GO" "$MATCHED_GO"

echo ">>> Job finished: $(date)"
