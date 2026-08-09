"""
GO Enrichment Analysis using gseapy (offline/local gene sets)
Input:  year_sig_gene_list.txt
Output: GO_ENRICHMENT/ directory
"""

import gseapy as gp
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os

GENE_FILE = "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/CLIMATE_ANALYSIS/year_sig_gene_list.txt"
OUT_DIR   = "/arf/scratch/ssenkal/selin/sekanslar/all_years_26/GO_ENRICHMENT"
os.makedirs(OUT_DIR, exist_ok=True)

# Read gene list
with open(GENE_FILE) as f:
    genes = [line.strip() for line in f if line.strip()]
print(f"Input genes: {len(genes)}")

# Local gene sets for Fly
gene_sets = {
    "BP": "GO_Biological_Process_2018",
    "MF": "GO_Molecular_Function_2018",
    "CC": "GO_Cellular_Component_2018",
    "KEGG": "KEGG_2019",
}

all_results = []

for cat, gs in gene_sets.items():
    print(f"\n>>> Running enrichment: {gs}")
    try:
        enr = gp.enrichr(
            gene_list     = genes,
            gene_sets     = gs,
            organism      = "fly",
            outdir        = None,
            cutoff        = 0.05,
            no_plot       = True,
        )

        if enr.results is not None and len(enr.results) > 0:
            sig = enr.results[enr.results["Adjusted P-value"] < 0.05].copy()
            sig["Ontology"] = cat
            print(f"  Significant terms: {len(sig)}")

            if len(sig) > 0:
                all_results.append(sig)
                sig.to_csv(f"{OUT_DIR}/GO_{cat}_results.tsv", sep="\t", index=False)

                # Dot plot
                top20 = sig.head(20).copy()
                top20["gene_count"] = top20["Overlap"].apply(lambda x: int(x.split("/")[0]))
                top20["-log10(padj)"] = -np.log10(top20["Adjusted P-value"].clip(lower=1e-300))

                fig, ax = plt.subplots(figsize=(11, max(6, len(top20)*0.45)))
                scatter = ax.scatter(
                    top20["-log10(padj)"],
                    range(len(top20)),
                    s=top20["gene_count"] * 25,
                    c=top20["Adjusted P-value"],
                    cmap="RdYlGn_r",
                    alpha=0.85,
                    edgecolors="grey",
                    linewidth=0.5
                )
                ax.set_yticks(range(len(top20)))
                ax.set_yticklabels(top20["Term"], fontsize=9)
                ax.set_xlabel("-log10(Adjusted P-value)", fontsize=11)
                ax.set_title(f"GO {cat} Enrichment\n(Top {len(top20)} significant terms, q<0.05)",
                             fontsize=12, fontweight="bold")
                cbar = plt.colorbar(scatter, ax=ax)
                cbar.set_label("Adjusted P-value", fontsize=9)
                plt.tight_layout()
                plt.savefig(f"{OUT_DIR}/dotplot_{cat}.pdf", bbox_inches="tight")
                plt.savefig(f"{OUT_DIR}/dotplot_{cat}.png", dpi=150, bbox_inches="tight")
                plt.close()

                # Bar plot
                fig, ax = plt.subplots(figsize=(11, max(6, len(top20)*0.45)))
                colors = plt.cm.RdYlGn_r(
                    (top20["Adjusted P-value"] - top20["Adjusted P-value"].min()) /
                    (top20["Adjusted P-value"].max() - top20["Adjusted P-value"].min() + 1e-10)
                )
                bars = ax.barh(range(len(top20)), top20["-log10(padj)"],
                               color=colors, edgecolor="grey", linewidth=0.5)
                ax.set_yticks(range(len(top20)))
                ax.set_yticklabels(top20["Term"], fontsize=9)
                ax.set_xlabel("-log10(Adjusted P-value)", fontsize=11)
                ax.set_title(f"GO {cat} Enrichment\n(Top {len(top20)} significant terms, q<0.05)",
                             fontsize=12, fontweight="bold")
                ax.invert_yaxis()
                plt.tight_layout()
                plt.savefig(f"{OUT_DIR}/barplot_{cat}.pdf", bbox_inches="tight")
                plt.savefig(f"{OUT_DIR}/barplot_{cat}.png", dpi=150, bbox_inches="tight")
                plt.close()

                print(f"  Plots saved: dotplot_{cat}.pdf, barplot_{cat}.pdf")
        else:
            print(f"  No significant terms")

    except Exception as e:
        print(f"  Error: {e}")

# Combined
if all_results:
    combined = pd.concat(all_results, ignore_index=True)
    combined = combined.sort_values("Adjusted P-value")
    combined.to_csv(f"{OUT_DIR}/GO_all_significant.tsv", sep="\t", index=False)
    print(f"\nTotal significant GO terms: {len(combined)}")

    # Summary bubble plot - all categories
    top10_each = combined.groupby("Ontology").head(10).copy()
    top10_each["gene_count"] = top10_each["Overlap"].apply(lambda x: int(x.split("/")[0]))
    top10_each["-log10(padj)"] = -np.log10(top10_each["Adjusted P-value"].clip(lower=1e-300))

    color_map = {"BP": "#E74C3C", "MF": "#3498DB", "CC": "#2ECC71", "KEGG": "#9B59B6"}
    fig, ax = plt.subplots(figsize=(13, max(8, len(top10_each)*0.4)))
    for _, row in top10_each.iterrows():
        ax.scatter(row["-log10(padj)"], row["Term"],
                   s=row["gene_count"] * 30,
                   c=color_map.get(row["Ontology"], "grey"),
                   alpha=0.8, edgecolors="grey", linewidth=0.5)
    from matplotlib.lines import Line2D
    legend_elements = [Line2D([0],[0], marker='o', color='w',
                               markerfacecolor=v, markersize=10, label=k)
                       for k, v in color_map.items() if k in top10_each["Ontology"].values]
    ax.legend(handles=legend_elements, title="Ontology", fontsize=9)
    ax.set_xlabel("-log10(Adjusted P-value)", fontsize=11)
    ax.set_title("GO Enrichment Summary — Year-Significant SNP Genes\n(Top 10 per category, q<0.05)",
                 fontsize=12, fontweight="bold")
    plt.tight_layout()
    plt.savefig(f"{OUT_DIR}/GO_summary_bubble.pdf", bbox_inches="tight")
    plt.savefig(f"{OUT_DIR}/GO_summary_bubble.png", dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Summary bubble plot saved.")
else:
    print("\nNo significant terms found.")

print(f"\n>>> Done: all results in {OUT_DIR}")
