import argparse, json, os
from datetime import datetime

def read_text(path, default=""):
    try:
        with open(path) as f:
            return f.read().strip()
    except:
        return default

def parse_quast_tsv(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                out[parts[0]] = parts[1]
    return out

def parse_kraken_species(report_path, topn=50):
    species = []
    if not os.path.exists(report_path):
        return species
    with open(report_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 6 and parts[3].strip() == "S":
                try:
                    species.append({
                        "name": parts[5].strip(),
                        "reads": int(parts[2].strip()),
                        "percentage": float(parts[0].strip())
                    })
                except:
                    pass
    species.sort(key=lambda x: x["percentage"], reverse=True)
    return species[:topn]

def load_json(path):
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path) as f:
            return json.load(f)
    except:
        return None

def is_skipped_file(path):
    txt = read_text(path, "")
    return ("SKIPPED" in txt) or (txt.strip() == "")

def write_text_report(path, report_obj):
    lines = []
    lines.append(f"Pipeline report (bioinformatician)")
    lines.append(f"Generated at: {report_obj.get('generated_at','')}")
    lines.append(f"Sample: {report_obj.get('sample_id','')}")
    lines.append(f"Platform: {report_obj.get('platform','')}")
    lines.append(f"Profile: {report_obj.get('profile','')}")
    lines.append("")
    lines.append("=== Coverage / AMR decision ===")
    lines.append(f"Mean coverage: {report_obj.get('mean_coverage',0)}")
    lines.append(f"AMR mode: {report_obj.get('amr_mode','unknown')}")
    lines.append("")

    lines.append("=== Taxonomy (top species from Kraken2 report) ===")
    top = report_obj.get("taxonomy", {}).get("kraken2_top_species", []) or []
    if not top:
        lines.append("No Kraken2 species parsed.")
    else:
        for s in top[:20]:
            lines.append(f"- {s.get('name','?')}: {s.get('percentage','?')}% ({s.get('reads','?')} reads)")
    lines.append("")

    lines.append("=== Centrifuge ===")
    cent = report_obj.get("taxonomy", {}).get("centrifuge", {}) or {}
    if not cent:
        lines.append("No centrifuge info.")
    else:
        lines.append(f"centrifuge_output: {cent.get('output_path','')}")
        lines.append(f"centrifuge_report: {cent.get('report_path','')}")
        lines.append(f"centrifuge_kreport: {cent.get('kreport_path','')}")
        lines.append(f"status: {cent.get('status','unknown')}")
    lines.append("")

    lines.append("=== Assembly (Quast summary TSV) ===")
    asm = report_obj.get("assembly", {}) or {}
    if not asm:
        lines.append("No assembly stats parsed.")
    else:
        for k in ["# contigs", "Total length", "Largest contig", "N50", "GC (%)"]:
            if k in asm:
                lines.append(f"{k}: {asm[k]}")
        shown = 0
        for k, v in asm.items():
            if k in {"# contigs", "Total length", "Largest contig", "N50", "GC (%)"}:
                continue
            lines.append(f"{k}: {v}")
            shown += 1
            if shown >= 15:
                break
    lines.append("")

    lines.append("=== AMR (contig-level) ===")
    amr = report_obj.get("amr", {}) or {}
    lines.append(f"ABRicate CARD path: {amr.get('abricate_card_path','')}")
    lines.append(f"AMRFinder path:     {amr.get('amrfinder_path','')}")
    lines.append("")

    bin_amr = report_obj.get("bin_amr", {}) or {}
    lines.append("=== Bin-level AMR / quality ===")
    if bin_amr:
        bq = bin_amr.get("bin_quality", {}) or {}
        if bq:
            lines.append(f"Bins with quality entries: {len(bq)}")
        else:
            lines.append("No bin_quality parsed.")
    else:
        lines.append("Not available (empty / missing).")
    lines.append("")

    warns = report_obj.get("warnings", []) or []
    if warns:
        lines.append("=== Warnings ===")
        for w in warns:
            lines.append(f"- {w}")
        lines.append("")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as f:
        f.write("\n".join(lines) + "\n")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--platform", required=True)
    ap.add_argument("--profile", required=True)
    ap.add_argument("--sample", required=True)

    ap.add_argument("--mean_cov", required=True)
    ap.add_argument("--amr_mode", required=True)

    ap.add_argument("--kraken_report", required=True)
    ap.add_argument("--bracken_report", required=True)

    ap.add_argument("--centrifuge_out", required=True)
    ap.add_argument("--centrifuge_report", required=True)
    ap.add_argument("--centrifuge_kreport", required=True)

    ap.add_argument("--quast_tsv", required=True)

    ap.add_argument("--abricate_card", required=True)
    ap.add_argument("--amrfinder", required=True)

    ap.add_argument("--bin_amr_json", required=True)

    ap.add_argument("--out", required=True)
    ap.add_argument("--out_txt", required=True)
    args = ap.parse_args()

    mean_cov = float(read_text(args.mean_cov, "0") or 0)
    amr_mode = read_text(args.amr_mode, "unknown")

    cent_status = "ok"
    if is_skipped_file(args.centrifuge_out) or is_skipped_file(args.centrifuge_report):
        cent_status = "skipped_or_missing_db"

    report = {
        "pipeline_version": "snakemake_v1",
        "generated_at": datetime.now().isoformat(),
        "sample_id": args.sample,
        "platform": args.platform,
        "profile": args.profile,
        "mean_coverage": mean_cov,
        "amr_mode": amr_mode,
        "taxonomy": {
            "kraken2_top_species": parse_kraken_species(args.kraken_report, 50),
            "bracken_report_path": args.bracken_report,
            "centrifuge": {
                "output_path": args.centrifuge_out,
                "report_path": args.centrifuge_report,
                "kreport_path": args.centrifuge_kreport,
                "status": cent_status
            }
        },
        "assembly": parse_quast_tsv(args.quast_tsv),
        "amr": {
            "abricate_card_path": args.abricate_card,
            "amrfinder_path": args.amrfinder
        },
        "warnings": []
    }

    report["bin_amr"] = load_json(args.bin_amr_json) or {}

    if amr_mode == "contig_low_nano":
        report["warnings"].append("LOW CONFIDENCE: Nanopore coverage <10x — assembly may be fragmented")
    if amr_mode == "reads_illumina":
        report["warnings"].append("LOW COVERAGE: Illumina <10x — KMA+ResFinder fallback may be used; contig AMR is backup.")

    if report.get("bin_amr") and report["bin_amr"].get("unbinned_amr_hits", 0) > 0:
        report["warnings"].append(
            f"UNBINNED: {report['bin_amr']['unbinned_amr_hits']} AMR hit(s) could not be assigned to any bin."
        )

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(report, f, indent=2)

    write_text_report(args.out_txt, report)

if __name__ == "__main__":
    main()
