import argparse
import json
import os
import sys
import shutil
import warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from sklearn.ensemble import RandomForestClassifier, IsolationForest
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.mixture import GaussianMixture
from sklearn.metrics import (ConfusionMatrixDisplay, RocCurveDisplay,
                             confusion_matrix, roc_auc_score)

def safe_float(val, default=0.0):
    try:
        s = str(val).strip()
        if '/' in s:
            num, denom = s.split('/', 1)
            denom = float(denom)
            return (float(num) / denom * 100.0) if denom != 0 else default
        return float(s) if s else default
    except (ValueError, TypeError):
        return default

# ─────────────────────────────────────────────
# GMM KMA to Host Linkage (DECLUTTERED GRAPH)
# ─────────────────────────────────────────────
def plot_kma_linkage(df_valid, kma_df, outpath):
    if df_valid.empty or kma_df.empty: return
    
    # 1. Filter organisms: Top 10 by depth + any organism that has a linked gene
    hosts_with_genes = set(kma_df["predicted_host"].unique())
    top_orgs = df_valid.sort_values("est_depth", ascending=False).head(10)["organism"].tolist()
    plot_orgs = list(set(top_orgs) | hosts_with_genes)
    
    df_plot = df_valid[df_valid["organism"].isin(plot_orgs)].sort_values("est_depth", ascending=False)
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    # Plot Organisms
    ax.scatter(df_plot["organism"], df_plot["est_depth"], c='lightblue', s=250, alpha=0.8, edgecolors='k', label="Host Organisms")
    
    # Plot KMA Genes with different colors/markers via a Legend
    top_kma = kma_df.sort_values("template_depth", ascending=False).head(10)
    genes = top_kma["gene"].values
    gene_depths = top_kma["template_depth"].values
    pred_hosts = top_kma["predicted_host"].values
    
    unique_genes = list(set(genes))
    cmap = plt.cm.get_cmap('tab10', max(1, len(unique_genes)))
    
    for i, g in enumerate(unique_genes):
        idx = np.where(genes == g)[0]
        # Plot the star marker
        ax.scatter(pred_hosts[idx], gene_depths[idx], color=cmap(i), marker='*', s=350, edgecolors='k', label=g, zorder=3)
        
        # Draw link lines
        for j in idx:
            h = pred_hosts[j]
            if h in df_plot["organism"].values:
                host_depth = df_plot[df_plot["organism"] == h]["est_depth"].values[0]
                ax.plot([h, h], [host_depth, gene_depths[j]], color=cmap(i), linestyle='--', alpha=0.7, zorder=2)

    plt.xticks(rotation=45, ha='right', fontsize=9)
    plt.ylabel("Sequencing Depth (x)")
    plt.title("Statistical Linkage: Pooling Raw-Read AMR to Host Genomes")
    
    # Move legend cleanly outside the plot
    ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", title="Detected AMR Genes", fontsize=9)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()

def link_kma_to_host(df: pd.DataFrame, kma_df: pd.DataFrame, outdir: str) -> pd.DataFrame:
    if kma_df.empty or df.empty:
        kma_df["predicted_host"] = "Unknown"
        return kma_df
        
    df_valid = df[df["rel_abundance"] > 0].copy()
    if df_valid.empty:
        kma_df["predicted_host"] = "Unknown"
        return kma_df
        
    mean_cov = df["mean_sample_coverage"].max() if "mean_sample_coverage" in df.columns else 10.0
    df_valid["est_depth"] = (df_valid["rel_abundance"] / 100.0) * mean_cov * 15.0 
    org_depths = df_valid["est_depth"].values.reshape(-1, 1)
    
    n_comp = min(3, len(org_depths))
    gmm = GaussianMixture(n_components=n_comp, random_state=42)
    gmm.fit(org_depths)
    
    predicted_hosts = []
    for _, row in kma_df.iterrows():
        closest_org = df_valid.iloc[(df_valid["est_depth"] - row["template_depth"]).abs().argsort()[:1]]["organism"].values[0]
        predicted_hosts.append(closest_org)
        
    kma_df["predicted_host"] = predicted_hosts
    plot_kma_linkage(df_valid, kma_df, os.path.join(outdir, "kma_host_linkage.png"))
    return kma_df

# ─────────────────────────────────────────────
# Feature extraction (FIXED KMA SEARCH)
# ─────────────────────────────────────────────
def extract_features(report: dict) -> pd.DataFrame:
    rows = []
    tax = report.get("taxonomy", {})
    top_species = tax.get("kraken2_top_species", []) or []
    bin_amr_data = report.get("bin_amr", {}) or {}
    bin_amr_hits = bin_amr_data.get("bin_amr_hits", []) or []
    bin_quality  = bin_amr_data.get("bin_quality", {}) or {}
    bin_orgs     = bin_amr_data.get("bin_organisms", {}) or {}
    plasmids     = report.get("plasmids", {}) or {}
    mean_cov     = float(report.get("mean_coverage", 0) or 0)
    
    # FIXED: Check BOTH locations for KMA hits to prevent them from vanishing
    amr_block = report.get("amr", {}) or {}
    kma_hits = amr_block.get("kma", {}).get("hits", []) or amr_block.get("reads_fallback", {}).get("hits", [])

    percs = [s.get("percentage", 0) for s in top_species if s.get("percentage", 0) > 0]
    p = np.array(percs) / 100.0
    alpha_div = float(-np.sum(p * np.log(p))) if len(p) > 0 else 0.0

    kma_n          = len(kma_hits)
    kma_mean_id    = float(np.mean([h.get("identity", 0) for h in kma_hits])) if kma_hits else 0.0
    kma_mean_dep   = float(np.mean([h.get("depth", 0) for h in kma_hits])) if kma_hits else 0.0
    kma_mean_cov   = float(np.mean([h.get("coverage_pct", 0) for h in kma_hits])) if kma_hits else 0.0

    has_plasmid = 1 if plasmids.get("count", 0) > 0 else 0

    for sp in top_species:
        name  = sp.get("name", "unknown").strip()
        pct   = float(sp.get("percentage", 0) or 0)
        reads = int(sp.get("reads", 0) or 0)

        matched_bins = [bname for bname, borg in bin_orgs.items() if isinstance(borg, str) and name.lower() in borg.lower()]
        
        # Merge quality stats for all matched bins (average completeness/contamination)
        completeness = 0.0
        contamination = 100.0
        if matched_bins:
            comps = [float(bin_quality.get(b, {}).get("completeness", 0) or 0) for b in matched_bins]
            conts = [float(bin_quality.get(b, {}).get("contamination", 100) or 100) for b in matched_bins]
            completeness = float(np.mean(comps)) if comps else 0.0
            contamination = float(np.mean(conts)) if conts else 100.0

        org_hits = [h for h in bin_amr_hits if h.get("bin") in matched_bins] if matched_bins else []
        
        # EXTRACT PROPER STRING NAMES
        gene_list = sorted(list(set([str(h.get("gene", "")).strip() for h in org_hits if h.get("gene")])))

        rows.append({
            "organism":             name,
            "bin":                  ", ".join(matched_bins) if matched_bins else "unmatched",
            "rel_abundance":        pct,
            "abs_reads":            reads,
            "completeness":         completeness,
            "contamination":        contamination,
            "n_amr_genes":          len(org_hits),
            "amr_gene_names":       ", ".join(gene_list),
            "high_conf_amr":        sum(1 for h in org_hits if h.get("high_confidence")),
            "plasmid_present":      has_plasmid,
            "chromosomal_amr":      sum(1 for h in org_hits if h.get("bin") not in ("unbinned",)),
            "plasmid_amr":          sum(1 for h in org_hits if h.get("bin") == "unbinned"),
            "max_amr_identity":     float(max([safe_float(h.get("identity", 0)) for h in org_hits] or [0.0])),
            "mean_amr_coverage":    float(np.mean([safe_float(h.get("coverage", 0)) for h in org_hits] or [0.0])),
            "kma_n_hits":           kma_n,
            "kma_mean_identity":    kma_mean_id,
            "kma_mean_depth":       kma_mean_dep,
            "kma_mean_coverage_pct":kma_mean_cov,
            "sample_alpha_diversity": alpha_div,
            "mean_sample_coverage":   mean_cov,
        })

    if not rows: return pd.DataFrame()
    return pd.DataFrame(rows)

FEATURE_COLS = [
    "rel_abundance", "abs_reads", "completeness", "contamination",
    "n_amr_genes", "high_conf_amr", "plasmid_present", "chromosomal_amr", "plasmid_amr",
    "max_amr_identity", "mean_amr_coverage", "kma_n_hits", "kma_mean_identity", 
    "kma_mean_depth", "kma_mean_coverage_pct", "sample_alpha_diversity", "mean_sample_coverage",
]

def heuristic_labels(df: pd.DataFrame) -> np.ndarray:
    cond = ((df["rel_abundance"] >= 1.0) & (df["completeness"] >= 50.0) & 
            (df["n_amr_genes"] >= 1) & (df["contamination"] < 15.0))
    labels = cond.astype(int).values
    if labels.sum() == 0: labels[df["rel_abundance"].values.argmax()] = 1
    return labels

def score_kma_precision(report: dict) -> pd.DataFrame:
    bin_amr_hits = (report.get("bin_amr", {}) or {}).get("bin_amr_hits", []) or []
    
    # FIXED: Check BOTH locations for hits
    amr_block = report.get("amr", {}) or {}
    kma_hits = amr_block.get("kma", {}).get("hits", []) or amr_block.get("reads_fallback", {}).get("hits", [])

    if not kma_hits:
        return pd.DataFrame(columns=["gene", "template_identity", "template_depth", "template_coverage_pct", "kma_precision_score", "predicted_host"])

    contig_features, kma_features, kma_rows = [], [], []
    for h in bin_amr_hits:
        try: contig_features.append([float(h.get("identity", 0) or 0), 0.0, float(h.get("coverage", 0) or 0)])
        except: pass

    for h in kma_hits:
        try:
            kma_features.append([float(h.get("identity", 0) or 0), float(h.get("depth", 0) or 0), float(h.get("coverage_pct", 0) or 0)])
            kma_rows.append({"gene": h.get("gene", "unknown"), "template_identity": float(h.get("identity", 0) or 0),
                             "template_depth": float(h.get("depth", 0) or 0), "template_coverage_pct": float(h.get("coverage_pct", 0) or 0)})
        except: pass

    if not kma_features: return pd.DataFrame()

    X_kma = np.array(kma_features)
    if len(contig_features) >= 3:
        iso = IsolationForest(n_estimators=100, contamination=0.1, random_state=42).fit(np.array(contig_features))
        raw_scores = iso.score_samples(X_kma)
        s_min, s_max = raw_scores.min(), raw_scores.max()
        precision_scores = (raw_scores - s_min) / (s_max - s_min) if s_max > s_min else np.ones(len(raw_scores)) * 0.5
    else:
        precision_scores = np.clip((X_kma[:, 0] / 100.0) * (X_kma[:, 2] / 100.0), 0, 1)

    kma_df = pd.DataFrame(kma_rows)
    kma_df["kma_precision_score"] = precision_scores
    return kma_df

def train_and_predict(df: pd.DataFrame, labels: np.ndarray, model_type: str = "rf") -> tuple[np.ndarray, object, pd.Series]:
    X = df[FEATURE_COLS].fillna(0).values
    X_scaled = StandardScaler().fit_transform(X)

    if model_type == "xgboost":
        try:
            from xgboost import XGBClassifier
            model = XGBClassifier(n_estimators=100, max_depth=3, learning_rate=0.1, use_label_encoder=False, eval_metric="logloss", random_state=42, n_jobs=4)
        except ImportError:
            model_type = "rf"

    if model_type == "rf":
        model = RandomForestClassifier(n_estimators=500, max_depth=None, min_samples_leaf=1, class_weight="balanced", random_state=42, n_jobs=4)

    model.fit(X_scaled, labels)
    proba = model.predict_proba(X_scaled)[:, 1]
    importances = pd.Series(model.feature_importances_, index=FEATURE_COLS) if hasattr(model, "feature_importances_") else pd.Series(np.zeros(len(FEATURE_COLS)), index=FEATURE_COLS)

    return proba, model, importances

# Plotting helpers omitted for brevity but they are preserved from previous script exactly
PALETTE = plt.cm.RdYlGn
def _save(fig, path): fig.tight_layout(); fig.savefig(path, dpi=150, bbox_inches="tight"); plt.close(fig)
def plot_empty(outpath: str, msg: str = "Not enough data"):
    fig, ax = plt.subplots(figsize=(6, 4)); ax.text(0.5, 0.5, msg, ha='center', va='center', fontsize=12, color='gray'); ax.axis('off'); _save(fig, outpath)
FEATURE_LABELS = { "rel_abundance": "Relative Abundance (%)", "abs_reads": "Total Reads", "completeness": "Genome Completeness (%)", "contamination": "Genome Contamination (%)", "n_amr_genes": "Total AMR Genes", "high_conf_amr": "High Confidence AMR Genes", "plasmid_present": "Plasmid Detected", "chromosomal_amr": "Chromosomal AMR Genes", "plasmid_amr": "Plasmid AMR Genes", "max_amr_identity": "Max AMR Identity (%)", "mean_amr_coverage": "Mean AMR Coverage (%)", "kma_n_hits": "KMA Fallback: AMR Gene Count", "kma_mean_identity": "KMA Fallback: Mean Identity (%)", "kma_mean_depth": "KMA Fallback: Mean Depth", "kma_mean_coverage_pct": "KMA Fallback: Mean Coverage (%)", "sample_alpha_diversity": "Sample Diversity (Shannon H)", "mean_sample_coverage": "Overall Sample Coverage" }
def _get_lbl(f): return FEATURE_LABELS.get(f, f)
def plot_feature_importance(importances: pd.Series, outpath: str):
    top = importances.sort_values(ascending=False).head(15); fig, ax = plt.subplots(figsize=(9, 5)); colors = PALETTE(np.linspace(0.7, 0.2, len(top))); top.index = [_get_lbl(x) for x in top.index]; top.plot(kind="barh", ax=ax, color=colors, edgecolor="k"); ax.set_xlabel("Importance Score (Influence on Prediction)"); ax.set_title("What Drives the Pathogen Prediction?"); ax.invert_yaxis(); _save(fig, outpath)
def plot_pca(df: pd.DataFrame, scores: np.ndarray, outpath: str):
    X = StandardScaler().fit_transform(df[FEATURE_COLS].fillna(0).values); pca = PCA(n_components=2, random_state=42); coords = pca.fit_transform(X); fig, ax = plt.subplots(figsize=(8, 6)); abundances = df["rel_abundance"].fillna(0).values; sizes = 20 + (abundances / (abundances.max() + 1e-9)) * 400; sc = ax.scatter(coords[:, 0], coords[:, 1], c=scores, cmap="RdYlGn", s=sizes, vmin=0, vmax=1, edgecolors="k", linewidths=0.5, alpha=0.9); df_sorted = df.sort_values(by="rel_abundance", ascending=False); top_3_indices = df_sorted.index[:3]
    for i, name in enumerate(df["organism"]):
        if scores[i] > 0.3 or i in top_3_indices: ax.annotate(name[:25], (coords[i, 0], coords[i, 1]), fontsize=7, ha="left", va="bottom", xytext=(3, 3), textcoords="offset points")
    cbar = plt.colorbar(sc, ax=ax); cbar.set_label("Likelihood of being True Pathogen (0-1)", rotation=270, labelpad=15); ax.set_xlabel(f"Variance Profile 1 ({pca.explained_variance_ratio_[0]*100:.1f}%)"); ax.set_ylabel(f"Variance Profile 2 ({pca.explained_variance_ratio_[1]*100:.1f}%)"); ax.set_title("Organism Clustering: Pathogens vs. Normal Flora"); _save(fig, outpath)
def plot_amr_heatmap(df: pd.DataFrame, outpath: str):
    amr_cols = ["n_amr_genes", "high_conf_amr", "plasmid_present", "chromosomal_amr", "plasmid_amr", "max_amr_identity"]; sub = df[["organism"] + amr_cols].set_index("organism"); fig, ax = plt.subplots(figsize=(8, max(4, len(sub) * 0.4))); im = ax.imshow(sub.values.astype(float), aspect="equal", cmap="YlOrRd"); plt.colorbar(im, ax=ax, label="Value", fraction=0.046, pad=0.04); ax.set_xticks(range(len(amr_cols))); ax.set_xticklabels([_get_lbl(c) for c in amr_cols], rotation=45, ha="right", fontsize=9); ax.set_yticks(range(len(sub))); ax.set_yticklabels([idx[:35] for idx in sub.index], fontsize=8); ax.set_title("Organism AMR Characteristics"); _save(fig, outpath)
def plot_kma_scatter(kma_df: pd.DataFrame, outpath: str):
    if kma_df.empty: plot_empty(outpath, "No KMA hits to score"); return
    fig, ax = plt.subplots(figsize=(10, 6)); genes = kma_df["gene"].unique(); colors = cm.get_cmap("tab20", max(1, len(genes)))
    for i, g in enumerate(genes):
        subset = kma_df[kma_df["gene"] == g]; ax.scatter(subset["template_identity"], subset["template_depth"], color=colors(i), s=80, edgecolors="k", label=g[:20], alpha=0.8)
    ax.set_xlabel("Sequence Match / Identity (%)"); ax.set_ylabel("Read Depth"); ax.set_title("AMR Gene Detection Reliability (KMA Fallback)")
    if len(genes) > 0: ax.legend(bbox_to_anchor=(1.04, 1), loc="upper left", title="Detected Genes", fontsize=8, ncol=1 + (len(genes)//20))
    _save(fig, outpath)
def plot_predictions_bar(df: pd.DataFrame, scores: np.ndarray, outpath: str):
    order = np.argsort(scores)[::-1]; names = [df["organism"].iloc[i][:35] for i in order]; vals  = scores[order]; fig, ax = plt.subplots(figsize=(10, max(4, len(names) * 0.45))); ax.barh(names, vals, color=[PALETTE(v) for v in vals], edgecolor="k", linewidth=0.5); ax.axvline(0.5, color="gray", linestyle="--", linewidth=1, label="Threshold (0.5)"); ax.set_xlim(0, 1); ax.set_xlabel("Pathogen Probability Score"); ax.set_title("Likelihood of Being the Causative Pathogen"); ax.legend(); ax.invert_yaxis(); _save(fig, outpath)
def plot_confusion(labels, preds, outpath: str):
    fig, ax = plt.subplots(figsize=(5, 4)); ConfusionMatrixDisplay(confusion_matrix(labels, preds), display_labels=["Background", "Pathogen"]).plot(ax=ax); ax.set_title("Confusion Matrix (ground truth labels)"); _save(fig, outpath)
def plot_roc(labels, scores, outpath: str):
    try: fig, ax = plt.subplots(figsize=(6, 5)); RocCurveDisplay.from_predictions(labels, scores, ax=ax); ax.set_title(f"ROC Curve (AUC = {roc_auc_score(labels, scores):.3f})"); _save(fig, outpath)
    except Exception: pass
def plot_taxonomy_concordance(report: dict, outpath: str):
    k_dict = {s["name"]: s.get("percentage", 0) for s in report.get("taxonomy", {}).get("kraken2_top_species", [])}
    c_dict = {s["name"]: s.get("abundance", 0) * 100 for s in report.get("taxonomy", {}).get("centrifuge_top_species", [])}
    sorted_spp = sorted(list(set(k_dict.keys()).union(set(c_dict.keys()))), key=lambda x: max(k_dict.get(x,0), c_dict.get(x,0)), reverse=True)
    fig, ax = plt.subplots(figsize=(8, 8)); cmap = cm.get_cmap("tab10", 10)
    for i, sp in enumerate(sorted_spp):
        kv, cv = k_dict.get(sp, 0.0), c_dict.get(sp, 0.0); color = cmap(i) if i < 10 and (kv > 1.0 or cv > 1.0) else "grey"
        ax.scatter([kv], [cv], c=[color], s=60, alpha=0.8 if color != "grey" else 0.4, edgecolors="k", label=sp[:30] if color != "grey" else None)
    max_val = max(max(k_dict.values()) if k_dict else 0, max(c_dict.values()) if c_dict else 0)
    ax.plot([0, max_val], [0, max_val], "k--", alpha=0.3, label="Perfect Concordance"); ax.set_xlabel("Kraken2 Abundance (%)"); ax.set_ylabel("Centrifuge Abundance (%)"); ax.set_title("Taxonomy Concordance (Kraken2 vs Centrifuge)")
    handles, labels = ax.get_legend_handles_labels(); by_label = dict(zip(labels, handles))
    if by_label: ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.04, 1), loc="upper left", title="Top Organisms", fontsize=8)
    _save(fig, outpath)
def plot_amr_concordance(report: dict, outpath: str):
    amr = report.get("amr", {}).get("contigs", {}).get("summary", {}); fig, ax = plt.subplots(figsize=(8, 8)); cmap = cm.get_cmap("tab20", 20); color_idx = 0
    for hit in amr.get("common", []): ax.scatter([hit.get("ab_identity", 0.0)], [hit.get("amr_identity", 0.0)], c=[cmap(color_idx % 20)], s=80, alpha=0.8, edgecolors="k", label=hit.get("gene", "unknown")[:20]); color_idx += 1
    for hit in amr.get("disjoint", []):
        src = hit.get("source", ""); ax.scatter([hit.get("identity", 0.0) if src == "ABRicate" else 0.0], [hit.get("identity", 0.0) if src == "AMRFinder" else 0.0], c=[cmap(color_idx % 20)], s=80, alpha=0.8, edgecolors="k", label=hit.get("gene", "unknown")[:20]); color_idx += 1
    ax.plot([0, 100], [0, 100], "k--", alpha=0.3, label="Perfect Concordance"); ax.set_xlabel("ABRicate Identity (%)"); ax.set_ylabel("AMRFinder Identity (%)"); ax.set_title("AMR Gene Concordance (Identity %)"); ax.set_xlim(-5, 105); ax.set_ylim(-5, 105)
    handles, labels = ax.get_legend_handles_labels(); by_label = dict(zip(labels, handles))
    if by_label: ax.legend(by_label.values(), by_label.keys(), bbox_to_anchor=(1.04, 1), loc="upper left", title="Detected Genes", fontsize=8, ncol=max(1, len(by_label) // 20))
    _save(fig, outpath)

# ─────────────────────────────────────────────
# Main
# ─────────────────────────────────────────────
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report_json", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--labels", default="")
    ap.add_argument("--model", default="rf", choices=["rf", "xgboost"])
    ap.add_argument("--pathogen_threshold", type=float, default=0.5)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    with open(args.report_json) as f: report = json.load(f)

    df = extract_features(report)
    if df.empty: sys.exit(1)

    kma_df = score_kma_precision(report)
    kma_df = link_kma_to_host(df, kma_df, args.outdir)

    if args.labels and os.path.exists(args.labels):
        gt_map = dict(zip(pd.read_csv(args.labels)["organism"], pd.read_csv(args.labels)["label"].astype(int)))
        labels = np.array([gt_map.get(o, 0) for o in df["organism"]])
        has_gt = True
    else:
        labels = heuristic_labels(df)
        has_gt = False

    scores, model, importances = train_and_predict(df, labels, args.model)

    out_df = df[["organism", "bin", "rel_abundance", "n_amr_genes", "amr_gene_names", "high_conf_amr", "plasmid_present", "completeness", "contamination"]].copy()
    out_df["pathogen_score"] = scores
    out_df["pathogen_flag"]  = (scores >= args.pathogen_threshold).astype(int)
    out_df["heuristic_label"] = labels
    out_df = out_df.sort_values("pathogen_score", ascending=False)
    
    # ── FIXED: Combine Contig Genes and ML-Linked Raw-Read Genes ──
    out_df["linked_amr_genes"] = ""
    out_df["kma_hit_count"] = 0
    if not kma_df.empty:
        linked_grouped = kma_df.groupby("predicted_host")["gene"].apply(lambda x: ", ".join(set(x))).to_dict()
        kma_counts = kma_df.groupby("predicted_host").size().to_dict()
        out_df["linked_amr_genes"] = out_df["organism"].map(linked_grouped).fillna("")
        out_df["kma_hit_count"] = out_df["organism"].map(kma_counts).fillna(0)
        
    def combine_genes(row):
        genes = set()
        if pd.notnull(row.get("amr_gene_names")) and row["amr_gene_names"]:
            genes.update([g.strip() for g in str(row["amr_gene_names"]).split(",") if g.strip()])
        if pd.notnull(row.get("linked_amr_genes")) and row["linked_amr_genes"]:
            genes.update([g.strip() for g in str(row["linked_amr_genes"]).split(",") if g.strip()])
        return ", ".join(sorted(list(genes)))

    out_df["final_amr_names"] = out_df.apply(combine_genes, axis=1)
    out_df["total_amr_hits"] = out_df["n_amr_genes"] + out_df["kma_hit_count"]

    pred_path = os.path.join(args.outdir, "predictions.csv")
    out_df.drop(columns=["final_amr_names", "linked_amr_genes", "kma_hit_count", "total_amr_hits"]).to_csv(pred_path, index=False)

    if not kma_df.empty:
        kma_df.to_csv(os.path.join(args.outdir, "kma_precision_scores.csv"), index=False)

    plot_feature_importance(importances, os.path.join(args.outdir, "feature_importance.png"))
    if len(df) >= 2: plot_pca(df, scores, os.path.join(args.outdir, "pca_clusters.png"))
    else: plot_empty(os.path.join(args.outdir, "pca_clusters.png"), "Not enough organisms for PCA")
        
    plot_amr_heatmap(df, os.path.join(args.outdir, "amr_heatmap.png"))
    plot_kma_scatter(kma_df, os.path.join(args.outdir, "kma_precision_scatter.png"))
    plot_predictions_bar(df, scores, os.path.join(args.outdir, "pathogen_scores_bar.png"))
    plot_taxonomy_concordance(report, os.path.join(args.outdir, "taxonomy_concordance.png"))
    plot_amr_concordance(report, os.path.join(args.outdir, "amr_concordance.png"))

    if has_gt:
        plot_confusion(labels, (scores >= args.pathogen_threshold).astype(int), os.path.join(args.outdir, "confusion_matrix.png"))
        plot_roc(labels, scores, os.path.join(args.outdir, "roc_curve.png"))

    ml_summary = {
        "model_type": args.model,
        "n_organisms_scored": len(out_df),
        "n_flagged_pathogens": int(out_df["pathogen_flag"].sum()),
        
        # FIXED: Pass the fully merged gene string and accurate total hit counts!
        "top_pathogens": out_df[out_df["pathogen_flag"] == 1][
            ["organism", "pathogen_score", "total_amr_hits", "final_amr_names"]
        ].rename(columns={"total_amr_hits": "n_amr_genes", "final_amr_names": "amr_gene_names"}).to_dict(orient="records"),
        
        "kma_n_scored": len(kma_df),
        "kma_mean_precision": float(kma_df["kma_precision_score"].mean()) if not kma_df.empty else None,
        "outputs": {
            "predictions_csv": pred_path,
            "kma_host_linkage": os.path.join(args.outdir, "kma_host_linkage.png")
        }
    }
    
    # NEW: Copy plasmid graph into ml_results for easy viewing
    plasmid_graph = report.get("plasmid_graph_path")
    if plasmid_graph and os.path.exists(plasmid_graph):
        dest = os.path.join(args.outdir, "plasmid_kmer_similarity.png")
        try:
            shutil.copy(plasmid_graph, dest)
            ml_summary["outputs"]["plasmid_graph"] = dest
        except: pass

    report["ml_results"] = ml_summary
    
    # Ensure KMA hits are saved back to the JSON properly
    if not kma_df.empty:
        report.setdefault("amr", {}).setdefault("kma", {})["linked_predictions"] = kma_df[["gene", "predicted_host", "kma_precision_score"]].to_dict('records')
        
    with open(args.report_json, "w") as f: json.dump(report, f, indent=2)

if __name__ == "__main__": main()
