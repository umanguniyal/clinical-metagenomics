import argparse
import json
import os
from datetime import datetime

TOPN_GENES = 20
TOPN_CLASSES = 20

def read_text(path, default=""):
    try:
        with open(path) as f:
            return f.read().strip()
    except:
        return default

def load_lines(path):
    try:
        with open(path) as f:
            return [ln.rstrip("\n") for ln in f]
    except:
        return []

def is_skipped_file(path):
    txt = read_text(path, "")
    return ("SKIPPED" in txt) or (txt.strip() == "")

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

def parse_abricate_tsv(path):
    hits = []
    if not os.path.exists(path) or is_skipped_file(path):
        return hits

    lines = load_lines(path)
    if not lines:
        return hits

    header = None
    for ln in lines:
        if ln.startswith("#") and "\t" in ln:
            header = ln.lstrip("#").strip().split("\t")
            break
        if ln.lower().startswith("file\t"):
            header = ln.strip().split("\t")
            break
    if header is None:
        header = lines[0].strip().split("\t")

    col = {name.strip().upper(): i for i, name in enumerate(header)}

    def get(fields, key, default=""):
        idx = col.get(key)
        if idx is None or idx >= len(fields):
            return default
        return fields[idx].strip()

    for ln in lines:
        if not ln or ln.startswith("#") or ln.lower().startswith("file\t"):
            continue
        fields = ln.split("\t")
        gene = get(fields, "GENE", "")
        if not gene and len(fields) > 5:
            gene = fields[5].strip()
        resistance = get(fields, "RESISTANCE", "")
        hits.append({"gene": gene or "unknown", "resistance": resistance or ""})

    return hits

def parse_amrfinder_tsv(path):
    hits = []
    if not os.path.exists(path) or is_skipped_file(path):
        return hits

    lines = load_lines(path)
    if not lines:
        return hits

    header = lines[0].split("\t")
    col = {name.strip().lower(): i for i, name in enumerate(header)}

    def get(fields, key, default=""):
        idx = col.get(key)
        if idx is None or idx >= len(fields):
            return default
        return fields[idx].strip()

    for ln in lines[1:]:
        if not ln.strip():
            continue
        fields = ln.split("\t")
        gene = get(fields, "gene symbol", "") or get(fields, "gene", "")
        if not gene and len(fields) >= 7:
            gene = fields[6].strip()
        drug_class = get(fields, "class", "") or get(fields, "subclass", "")
        hits.append({"gene": gene or "unknown", "class": drug_class or ""})

    return hits

def summarize_hits(hits, gene_key="gene", class_key=None, topn_genes=TOPN_GENES, topn_classes=TOPN_CLASSES):
    gene_counts = {}
    class_counts = {}
    for h in hits:
        g = (h.get(gene_key) or "unknown").strip()
        gene_counts[g] = gene_counts.get(g, 0) + 1
        if class_key:
            c = (h.get(class_key) or "").strip()
            if c:
                class_counts[c] = class_counts.get(c, 0) + 1

    top_genes = sorted(gene_counts.items(), key=lambda x: (-x[1], x[0]))[:topn_genes]
    top_classes = sorted(class_counts.items(), key=lambda x: (-x[1], x[0]))[:topn_classes]

    return {
        "total_hits": len(hits),
        "unique_genes": len(gene_counts),
        "top_genes": [{"gene": g, "count": n} for g, n in top_genes],
        "resistance_classes": [{"class": c, "count": n} for c, n in top_classes],
    }

def kma_res_summary(path, max_lines=30):
    if not path or not os.path.exists(path) or is_skipped_file(path):
        return {"ran": False, "preview": []}
    lines = load_lines(path)
    filtered = [ln for ln in lines if ln.strip() and "SKIPPED" not in ln]
    return {"ran": True, "preview": filtered[:max_lines]}

def write_text_report(path, report_obj):
    lines = []
    lines.append(f"Pipeline report (clinician)")
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

    lines.append("=== AMR summary ===")
    amr = report_obj.get("amr", {}) or {}
    cont = (amr.get("contigs") or {})
    reads = (amr.get("reads_fallback") or {})

    lines.append("-- Contig-based AMR (ABRicate + AMRFinder) --")
    lines.append(f"ABRicate CARD: {cont.get('abricate_card_path','')}")
    lines.append(f"AMRFinder:     {cont.get('amrfinder_path','')}")
    cont_summary = cont.get("summary") or {}
    if cont_summary:
        lines.append(f"Total contig AMR hits: {cont_summary.get('total_hits',0)}")
        lines.append(f"Unique genes: {cont_summary.get('unique_genes',0)}")
        tg = cont_summary.get("top_genes", []) or []
        if tg:
            lines.append("Top genes:")
            for x in tg[:TOPN_GENES]:
                lines.append(f"  - {x['gene']} (n={x['count']})")
        rc = cont_summary.get("resistance_classes", []) or []
        if rc:
            lines.append("Resistance classes:")
            for x in rc[:TOPN_CLASSES]:
                lines.append(f"  - {x['class']} (n={x['count']})")
        else:
            lines.append("Resistance classes: N/A")
    else:
        lines.append("Contig AMR summary: SKIPPED or not available.")
    lines.append("")

    lines.append("-- Read-based AMR fallback (KMA/ResFinder) --")
    if reads.get("kma_res_path"):
        lines.append(f"KMA .res: {reads.get('kma_res_path')}")
        prev = reads.get("kma_res_preview", []) or []
        if prev:
            lines.append("KMA preview:")
            for ln in prev:
                lines.append(f"  {ln}")
    else:
        lines.append("Read-based AMR: not run / SKIPPED.")
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

    ap.add_argument("--kma_res", default="")
    ap.add_argument("--kma_tsv", default="")

    ap.add_argument("--out", required=True)
    ap.add_argument("--out_txt", required=True)

    args = ap.parse_args()

    mean_cov = float(read_text(args.mean_cov, "0") or 0)
    amr_mode = read_text(args.amr_mode, "unknown")

    top_species = parse_kraken_species(args.kraken_report, 50)

    cent_status = "ok"
    if is_skipped_file(args.centrifuge_out) or is_skipped_file(args.centrifuge_report):
        cent_status = "skipped_or_missing_db"

    ab_hits = parse_abricate_tsv(args.abricate_card)
    amrf_hits = parse_amrfinder_tsv(args.amrfinder)

    combined_hits = [{"gene": h["gene"], "resistance": h.get("resistance","")} for h in ab_hits]
    combined_hits += [{"gene": h["gene"], "resistance": h.get("class","")} for h in amrf_hits]
    contig_summary = summarize_hits(combined_hits, gene_key="gene", class_key="resistance")

    kma_prev = kma_res_summary(args.kma_res).get("preview", [])
    kma_ran = bool(kma_prev) or (args.kma_res and os.path.exists(args.kma_res) and not is_skipped_file(args.kma_res))

    report = {
        "pipeline_version": "snakemake_v1",
        "generated_at": datetime.now().isoformat(),
        "sample_id": args.sample,
        "platform": args.platform,
        "profile": args.profile,
        "mean_coverage": mean_cov,
        "amr_mode": amr_mode,
        "taxonomy": {
            "kraken2_top_species": top_species,
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
            "contigs": {
                "abricate_card_path": args.abricate_card,
                "amrfinder_path": args.amrfinder,
                "summary": contig_summary
            },
            "reads_fallback": {}
        },
        "warnings": []
    }

    if kma_ran:
        report["amr"]["reads_fallback"] = {
            "kma_res_path": args.kma_res,
            "kma_tsv_path": args.kma_tsv,
            "kma_res_preview": kma_prev[:30]
        }

    if amr_mode == "reads_illumina":
        report["warnings"].append("LOW COVERAGE: Illumina < threshold — using read-based AMR (KMA/ResFinder).")
    if amr_mode == "contig_low_nano":
        report["warnings"].append("LOW CONFIDENCE: Nanopore coverage < threshold — contig AMR still reported (no read fallback).")

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        json.dump(report, f, indent=2)

    write_text_report(args.out_txt, report)

if __name__ == "__main__":
    main()
