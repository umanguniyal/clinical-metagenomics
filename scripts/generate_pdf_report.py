"""
generate_pdf_report.py

Reads a JSON report (clinician or bioinformatician) and renders a
publication-quality PDF using Jinja2 + WeasyPrint.

Usage:
python generate_pdf_report.py --report_json report.json --out report.pdf
"""
import argparse
import json
import os
import sys
from datetime import datetime

from jinja2 import Environment, FileSystemLoader
from weasyprint import HTML

TEMPLATE_DIR = os.path.join(os.path.dirname(__file__), "templates")


def load_report(path):
    with open(path) as f:
        return json.load(f)


def prepare_context(report):
    """Flatten the JSON into simple template variables."""
    ctx = {
        "sample_id":        report.get("sample_id", "unknown"),
        "generated_at":     report.get("generated_at", datetime.now().isoformat()),
        "pipeline_version": report.get("pipeline_version", "snakemake_v1"),
        "platform":         report.get("platform", "unknown"),
        "profile":          report.get("profile", "unknown"),
        "mean_coverage":    report.get("mean_coverage", 0),
        "amr_mode":         report.get("amr_mode", "unknown"),
        "assembly":         report.get("assembly", {}),
        "host_removal":     report.get("host_removal", {}),
        "warnings":         report.get("warnings", []),
    }

    # Taxonomy
    taxonomy = report.get("taxonomy", {})
    k_spp = taxonomy.get("kraken2_top_species", [])
    c_spp = taxonomy.get("centrifuge_top_species", [])
    
    # Merge taxonomy
    k_dict = {s["name"]: s.get("percentage", 0) for s in k_spp}
    c_dict = {s["name"]: s.get("abundance", 0) * 100 for s in c_spp} # Convert fraction to pct
    
    all_spp = set(k_dict.keys()).union(set(c_dict.keys()))
    merged_tax = []
    for sp in all_spp:
        merged_tax.append({
            "name": sp,
            "k_pct": k_dict.get(sp),
            "c_pct": c_dict.get(sp)
        })
    
    # Sort: Common first (by max abundance), then disjoint (by max abundance)
    merged_tax.sort(key=lambda x: (
        0 if (x["k_pct"] is not None and x["c_pct"] is not None) else 1,
        -max(x["k_pct"] or 0, x["c_pct"] or 0)
    ))
    ctx["merged_taxonomy"] = merged_tax[:40]

    # AMR
    amr = report.get("amr", {})
    contigs = amr.get("contigs", {})
    summary = contigs.get("summary", {})
    
    if "common" in summary or "disjoint" in summary:
        ctx["merged_amr"] = {
            "common": summary.get("common", []),
            "disjoint": summary.get("disjoint", [])
        }
    else:
        ctx["merged_amr"] = {"common": [], "disjoint": []}

    # --- KMA Status Information ---
    # Check both possible locations for KMA data
    kma_info = {}
    if amr.get("kma"):
        kma_info = amr["kma"]
    elif amr.get("reads_fallback", {}).get("kma_res_path"):
        kma_info = amr["reads_fallback"]
    
    ctx["kma_status"] = kma_info.get("status", "unknown")
    ctx["kma_hits"] = kma_info.get("hits", [])
    
    # --- ML Data Injections ---
    ctx["ml_results"] = report.get("ml_results", {})
    # Only keep ML-based KMA host linkages with a non-zero precision score
    _raw_linked = amr.get("kma", {}).get("linked_predictions", []) or []
    ctx["linked_predictions"] = [
        p for p in _raw_linked
        if float(p.get("kma_precision_score", 0) or 0) > 0
    ]

    # Plasmids
    ctx["plasmids"] = report.get("plasmids", {})

    # Binned AMR
    ctx["bin_amr"] = report.get("bin_amr", {})

    return ctx


def main():
    ap = argparse.ArgumentParser(description="Generate PDF report from JSON")
    ap.add_argument("--report_json", required=True, help="Path to the JSON report")
    ap.add_argument("--out", required=True, help="Output PDF path")
    args = ap.parse_args()

    report = load_report(args.report_json)
    ctx = prepare_context(report)

    env = Environment(loader=FileSystemLoader(TEMPLATE_DIR))
    template = env.get_template("report.html")
    html_str = template.render(**ctx)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    HTML(string=html_str).write_pdf(args.out)
    print(f"PDF written to: {args.out}", file=sys.stderr)


if __name__ == "__main__":
    main()