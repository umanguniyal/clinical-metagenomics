"""
aggregate_cohort_ml.py
======================
Generates cumulative multi-sample analysis graphs across the entire cohort.
"""

import argparse
import json
import os
import sys
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

FEATURE_COLS = [
    "rel_abundance", "abs_reads",
    "completeness", "contamination",
    "n_amr_genes", "high_conf_amr",
    "plasmid_present", "chromosomal_amr", "plasmid_amr",
    "max_amr_identity", "mean_amr_coverage",
    "kma_n_hits", "kma_mean_identity", "kma_mean_depth", "kma_mean_coverage_pct",
    "sample_alpha_diversity", "mean_sample_coverage",
]

def _save(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}", file=sys.stderr)

def extract_cohort_features(json_files):
    from ml_pathogen_amr import extract_features
    rows = []
    amr_associations = []
    concordances = []
    
    for path in json_files:
        if not os.path.exists(path):
            continue
        try:
            with open(path) as f:
                report = json.load(f)
        except Exception:
            continue
            
        sample_id = report.get("sample_id", os.path.basename(os.path.dirname(os.path.dirname(path))))
        
        # Get ML scores if they exist
        ml_results = report.get("ml_results", {})
        top_pathogens = {x["organism"]: x["pathogen_score"] for x in ml_results.get("top_pathogens", [])}
        
        # Extract features just like single sample
        df = extract_features(report)
        if df.empty:
            continue
            
        df["sample_id"] = sample_id
        df["pathogen_score"] = df["organism"].apply(lambda x: top_pathogens.get(x, 0.0))
        
        rows.append(df)
        
        # Extract Species-AMR association
        bin_amr = report.get("bin_amr", {}).get("bin_amr_hits", []) or []
        bin_orgs = report.get("bin_amr", {}).get("bin_organisms", {}) or {}
        
        for h in bin_amr:
            bname = h.get("bin")
            gene = h.get("gene")
            if not bname or not gene: continue
            
            org = "Unknown"
            if bname in bin_orgs and isinstance(bin_orgs[bname], str):
                org = bin_orgs[bname]
            
            amr_associations.append({
                "sample_id": sample_id,
                "organism": org,
                "gene": gene
            })
            
        tax = report.get("taxonomy", {})
        k_spp = tax.get("kraken2_top_species", [])
        c_spp = tax.get("centrifuge_top_species", [])
        k_dict = {s["name"]: s.get("percentage", 0) for s in k_spp}
        c_dict = {s["name"]: s.get("abundance", 0) * 100 for s in c_spp}
        tax_intersect = set(k_dict.keys()).intersection(set(c_dict.keys()))
        tax_union = set(k_dict.keys()).union(set(c_dict.keys()))
        tax_conc = len(tax_intersect) / len(tax_union) if tax_union else 0.0
        
        amr = report.get("amr", {}).get("contigs", {}).get("summary", {})
        amr_common = amr.get("common", [])
        amr_disjoint = amr.get("disjoint", [])
        amr_intersect = len(amr_common)
        amr_union = amr_intersect + len(amr_disjoint)
        amr_conc = amr_intersect / amr_union if amr_union > 0 else 1.0
        
        concordances.append({
            "sample_id": sample_id,
            "tax_concordance": tax_conc,
            "amr_concordance": amr_conc
        })
            
    if not rows:
        return pd.DataFrame(), pd.DataFrame(), pd.DataFrame()
        
    return pd.concat(rows, ignore_index=True), pd.DataFrame(amr_associations), pd.DataFrame(concordances)

def plot_cohort_pca(df, outpath):
    if len(df) < 2: return
    X = df[FEATURE_COLS].fillna(0).values
    X_s = StandardScaler().fit_transform(X)
    pca = PCA(n_components=2, random_state=42)
    coords = pca.fit_transform(X_s)

    fig, ax = plt.subplots(figsize=(10, 8))
    
    abundances = df["rel_abundance"].fillna(0).values
    sizes = 20 + (abundances / (abundances.max() + 1e-9)) * 400
    
    scores = df["pathogen_score"].values

    sc = ax.scatter(coords[:, 0], coords[:, 1],
                    c=scores, cmap="RdYlGn", s=sizes,
                    vmin=0, vmax=1, edgecolors="k", linewidths=0.5, alpha=0.8)
    
    cbar = plt.colorbar(sc, ax=ax)
    cbar.set_label("Pathogen Score (0=Flora, 1=Pathogen)", rotation=270, labelpad=15)
    
    ax.set_xlabel(f"Variance Profile 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"Variance Profile 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title(f"Cohort Landscape: Pathogens vs. Normal Flora (N={df['sample_id'].nunique()} samples)")
    _save(fig, outpath)

def plot_amr_burden(df, outpath):
    if df.empty: return
    
    # Aggregate AMR burden per sample
    grouped = df.groupby("sample_id")[["chromosomal_amr", "plasmid_amr"]].sum().reset_index()
    grouped = grouped.sort_values(by="plasmid_amr", ascending=False)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    ax.bar(grouped["sample_id"], grouped["chromosomal_amr"], label="Chromosomal AMR", color="#3498db", edgecolor="k")
    ax.bar(grouped["sample_id"], grouped["plasmid_amr"], bottom=grouped["chromosomal_amr"], label="Plasmid-Borne AMR", color="#e74c3c", edgecolor="k")
    
    ax.set_xlabel("Sample ID")
    ax.set_ylabel("Total AMR Genes Detected")
    ax.set_title("AMR Burden Tracker Across Cohort")
    ax.legend()
    plt.xticks(rotation=45, ha="right", fontsize=8)
    _save(fig, outpath)

def plot_species_amr_matrix(assoc_df, outpath):
    if assoc_df.empty: return
    
    # Create cross-tabulation
    ct = pd.crosstab(assoc_df["organism"], assoc_df["gene"])
    if ct.empty: return
    
    # Filter for top 30 organisms and top 30 genes by frequency to avoid massive plots
    top_orgs = ct.sum(axis=1).sort_values(ascending=False).head(30).index
    top_genes = ct.sum(axis=0).sort_values(ascending=False).head(30).index
    
    sub = ct.loc[top_orgs, top_genes]
    
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(sub.values, aspect="auto", cmap="Reds")
    plt.colorbar(im, ax=ax, label="Detection Frequency")
    
    ax.set_xticks(range(len(sub.columns)))
    ax.set_xticklabels(sub.columns, rotation=90, fontsize=8)
    ax.set_yticks(range(len(sub.index)))
    ax.set_yticklabels(sub.index, fontsize=8)
    ax.set_title("Species-AMR Association Matrix (Top 30)")
    
    _save(fig, outpath)

def plot_score_dist(df, outpath):
    if df.empty: return
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(df["pathogen_score"], bins=20, color="#2ecc71", edgecolor="k", alpha=0.7)
    ax.set_xlabel("Pathogen Probability Score")
    ax.set_ylabel("Frequency (Organisms)")
    ax.set_title("Distribution of Pathogen Scores Across Cohort")
    _save(fig, outpath)

def plot_top_pathogens(df, outpath):
    if df.empty: return
    
    # Filter for high confidence pathogens
    pathogens = df[df["pathogen_score"] >= 0.5]
    if pathogens.empty:
        # Write empty plot
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.text(0.5, 0.5, "No pathogens detected in cohort", ha='center', va='center', fontsize=12, color='gray')
        ax.axis('off')
        _save(fig, outpath)
        return
        
    counts = pathogens["organism"].value_counts().head(20)
    fig, ax = plt.subplots(figsize=(10, 6))
    
    counts.plot(kind="bar", color="#e74c3c", edgecolor="k", ax=ax)
    ax.set_xlabel("Organism")
    ax.set_ylabel("Number of Samples")
    ax.set_title("Top Flagged Pathogens Across Cohort")
    plt.xticks(rotation=45, ha="right", fontsize=9)
    _save(fig, outpath)

def plot_cohort_concordance(conc_df, outpath):
    if conc_df.empty: return
    
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.boxplot([conc_df["tax_concordance"] * 100, conc_df["amr_concordance"] * 100], labels=["Taxonomy (Kraken2 vs Centrifuge)", "AMR (ABRicate vs AMRFinder)"])
    
    ax.set_ylabel("Jaccard Concordance (%)")
    ax.set_title("Pipeline Internal Reliability Across Cohort")
    _save(fig, outpath)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reports", nargs="+", required=True, help="List of report_bioinformatician.json files")
    ap.add_argument("--outdir", required=True, help="Output directory for cohort plots")
    args = ap.parse_args()
    
    os.makedirs(args.outdir, exist_ok=True)
    
    print("Loading cohort data...", file=sys.stderr)
    sys.path.append(os.path.dirname(__file__)) # For importing ml_pathogen_amr
    
    df, assoc_df, conc_df = extract_cohort_features(args.reports)
    if df.empty:
        print("No valid data found in provided reports.", file=sys.stderr)
        sys.exit(0)
        
    print(f"Loaded {len(df)} organisms from {df['sample_id'].nunique()} samples.", file=sys.stderr)
    
    plot_cohort_pca(df, os.path.join(args.outdir, "cohort_pca_landscape.png"))
    plot_amr_burden(df, os.path.join(args.outdir, "cohort_amr_burden.png"))
    plot_species_amr_matrix(assoc_df, os.path.join(args.outdir, "cohort_species_amr_matrix.png"))
    plot_score_dist(df, os.path.join(args.outdir, "cohort_score_distribution.png"))
    plot_top_pathogens(df, os.path.join(args.outdir, "cohort_top_pathogens.png"))
    plot_cohort_concordance(conc_df, os.path.join(args.outdir, "cohort_concordance.png"))
    
    print("Done generating cumulative cohort graphs.", file=sys.stderr)

if __name__ == "__main__":
    main()
