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

def parse_fasta_lengths(path):
    lengths = []
    if not path or not os.path.exists(path):
        return lengths
    curr_len = 0
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if curr_len > 0:
                    lengths.append(curr_len)
                curr_len = 0
            else:
                curr_len += len(line.strip())
        if curr_len > 0:
            lengths.append(curr_len)
    return lengths

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

def parse_centrifuge_species(report_path, topn=50):
    species = []
    if not report_path or not os.path.exists(report_path):
        return species
    with open(report_path) as f:
        header = f.readline().rstrip('\n').split('\t')
        try:
            name_idx = header.index("name")
            rank_idx = header.index("taxRank")
            reads_idx = header.index("numReads")
            abund_idx = header.index("abundance")
        except ValueError:
            return species
            
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) > max(name_idx, rank_idx, reads_idx, abund_idx):
                if parts[rank_idx].strip() == "species":
                    try:
                        species.append({
                            "name": parts[name_idx].strip(),
                            "reads": int(parts[reads_idx].strip()),
                            "abundance": float(parts[abund_idx].strip())
                        })
                    except:
                        pass
    species.sort(key=lambda x: (x["reads"], x["abundance"]), reverse=True)
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
        contig = get(fields, "SEQUENCE", "")
        if not contig and len(fields) > 1:
            contig = fields[1].split()[0]
        
        # Clean up ABRicate gene names if needed
        gene = gene.split(".")[0] if "." in gene else gene
        identity = get(fields, "%IDENTITY", "0")
        coverage = get(fields, "%COVERAGE", "0")
        hits.append({"gene": gene or "unknown", "resistance": resistance or "", "contig": contig, "identity": identity, "coverage": coverage})
    return hits

def parse_kma_res(path):
    hits = []
    if not path or not os.path.exists(path) or is_skipped_file(path):
        return hits
    try:
        with open(path) as f:
            header = f.readline().rstrip('\n').lstrip('#').split('\t')
            header = [h.strip() for h in header]
            for line in f:
                if not line.strip() or line.startswith('SKIPPED'):
                    continue
                parts = line.rstrip('\n').split('\t')
                row = dict(zip(header, parts))
                try:
                    hits.append({
                        "gene":          row.get("# Template") or row.get("Template") or parts[0],
                        "identity":      float(row.get("Template_Identity", 0) or 0),
                        "coverage_pct":  float(row.get("Template_Coverage", 0) or 0),
                        "depth":         float(row.get("Depth", 0) or 0)
                    })
                except (ValueError, IndexError):
                    pass
    except Exception:
        pass
    hits.sort(key=lambda x: (x["depth"], x["identity"]), reverse=True)
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
        contig = get(fields, "contig id", "")
        if not contig and len(fields) > 1:
            contig = fields[1].split()[0]
        identity = get(fields, "% identity to reference sequence", "0")
        coverage = get(fields, "% coverage of reference sequence", "0")
        hits.append({"gene": gene or "unknown", "class": drug_class or "", "contig": contig, "identity": identity, "coverage": coverage})
    return hits

def summarize_hits(ab_hits, amrf_hits):
    ab_counts = {}
    for h in ab_hits:
        g = h.get("gene", "unknown")
        if g not in ab_counts:
            ab_counts[g] = {"count": 0, "res": h.get("resistance", ""), "identities": [], "coverages": []}
        ab_counts[g]["count"] += 1
        if not ab_counts[g]["res"]: ab_counts[g]["res"] = h.get("resistance", "")
        try:
            ab_counts[g]["identities"].append(float(h.get("identity", 0)))
            ab_counts[g]["coverages"].append(float(h.get("coverage", 0)))
        except:
            pass
            
    amr_counts = {}
    for h in amrf_hits:
        g = h.get("gene", "unknown")
        if g not in amr_counts:
            amr_counts[g] = {"count": 0, "res": h.get("class", ""), "identities": [], "coverages": []}
        amr_counts[g]["count"] += 1
        if not amr_counts[g]["res"]: amr_counts[g]["res"] = h.get("class", "")
        try:
            amr_counts[g]["identities"].append(float(h.get("identity", 0)))
            amr_counts[g]["coverages"].append(float(h.get("coverage", 0)))
        except:
            pass

    def median(lst):
        if not lst: return 0.0
        s = sorted(lst)
        n = len(s)
        return (s[n//2] + s[(n-1)//2]) / 2.0

    common_genes = set(ab_counts.keys()).intersection(amr_counts.keys())
    only_ab = set(ab_counts.keys()) - common_genes
    only_amr = set(amr_counts.keys()) - common_genes

    common = []
    for g in common_genes:
        ab_max_id = max(ab_counts[g]["identities"]) if ab_counts[g]["identities"] else 0.0
        amr_max_id = max(amr_counts[g]["identities"]) if amr_counts[g]["identities"] else 0.0
        ab_cov = median(ab_counts[g]["coverages"])
        amr_cov = median(amr_counts[g]["coverages"])
        common.append({
            "gene": g,
            "resistance": ab_counts[g]["res"] or amr_counts[g]["res"],
            "count": max(ab_counts[g]["count"], amr_counts[g]["count"]),
            "ab_identity": ab_max_id,
            "amr_identity": amr_max_id,
            "ab_coverage": ab_cov,
            "amr_coverage": amr_cov,
            "median_coverage": (ab_cov + amr_cov) / 2.0
        })
    common.sort(key=lambda x: -x["count"])

    disjoint = []
    for g in only_ab:
        ab_max_id = max(ab_counts[g]["identities"]) if ab_counts[g]["identities"] else 0.0
        ab_cov = median(ab_counts[g]["coverages"])
        disjoint.append({
            "gene": g, "resistance": ab_counts[g]["res"], 
            "count": ab_counts[g]["count"], "source": "ABRicate",
            "identity": ab_max_id,
            "coverage": ab_cov,
            "median_coverage": ab_cov
        })
    for g in only_amr:
        amr_max_id = max(amr_counts[g]["identities"]) if amr_counts[g]["identities"] else 0.0
        amr_cov = median(amr_counts[g]["coverages"])
        disjoint.append({
            "gene": g, "resistance": amr_counts[g]["res"], 
            "count": amr_counts[g]["count"], "source": "AMRFinder",
            "identity": amr_max_id,
            "coverage": amr_cov,
            "median_coverage": amr_cov
        })
    disjoint.sort(key=lambda x: -x["count"])

    return {
        "common": common,
        "disjoint": disjoint
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

    lines.append("=== Taxonomy (top species from Centrifuge) ===")
    top_cent = report_obj.get("taxonomy", {}).get("centrifuge_top_species", []) or []
    if not top_cent:
        lines.append("No Centrifuge species parsed.")
    else:
        for s in top_cent[:20]:
            lines.append(f"- {s.get('name','?')}: abundance={s.get('abundance',0):.5f} ({s.get('reads','?')} reads)")
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
        lines.append(f"Common AMR hits (both tools): {len(cont_summary.get('common', []))}")
        for x in cont_summary.get("common", []):
            lines.append(f"  - {x['gene']} [{x['resistance']}] (n={x['count']})")
        
        lines.append(f"Disjoint AMR hits (one tool): {len(cont_summary.get('disjoint', []))}")
        for x in cont_summary.get("disjoint", []):
            lines.append(f"  - {x['gene']} [{x['resistance']}] (n={x['count']}) - {x['source']}")
    else:
        lines.append("Contig AMR summary: SKIPPED or not available.")
    lines.append("")

    # REVERTED: No ML logic here. Just print the raw KMA preview.
    lines.append("-- Read-based AMR fallback (KMA/ResFinder) --")
    if reads.get("kma_res_path") or (amr.get("kma", {}).get("res_path")):
        lines.append(f"KMA output detected.")
        prev = reads.get("kma_res_preview", []) or amr.get("kma", {}).get("hits", [])
        if prev:
            lines.append("KMA preview:")
            for ln in prev[:15]:
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
    ap.add_argument("--plasmid_fasta", default="")
    ap.add_argument("--mobtyper_results", default="")
    ap.add_argument("--contig_report", default="")
    ap.add_argument("--bins_dir", default="")
    ap.add_argument("--host_stats", default="")

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

    contig_summary = summarize_hits(ab_hits, amrf_hits)

    kma_prev = kma_res_summary(args.kma_res).get("preview", [])
    kma_ran = bool(kma_prev) or (args.kma_res and os.path.exists(args.kma_res) and not is_skipped_file(args.kma_res))

    host_stats = {}
    if args.host_stats and os.path.exists(args.host_stats):
        try:
            with open(args.host_stats) as f:
                host_stats = json.load(f)
        except:
            pass

    report = {
        "pipeline_version": "snakemake_v1",
        "generated_at": datetime.now().isoformat(),
        "sample_id": args.sample,
        "platform": args.platform,
        "profile": args.profile,
        "mean_coverage": mean_cov,
        "amr_mode": amr_mode,
        "host_removal": host_stats,
        "taxonomy": {
            "kraken2_top_species": top_species,
            "centrifuge_top_species": parse_centrifuge_species(args.centrifuge_report, 50),
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

# KMA hits processing - always include for consistent reporting
    kma_hits = parse_kma_res(args.kma_res) if args.kma_res and os.path.exists(args.kma_res) else []
    kma_status = "completed" if kma_hits else ("skipped" if is_skipped_file(args.kma_res) or not os.path.exists(args.kma_res) else "no_hits")
    
    report["amr"]["reads_fallback"] = {
        "kma_res_path": args.kma_res,
        "kma_tsv_path": args.kma_tsv,
        "hits": kma_hits,
        "kma_res_preview": kma_prev[:30] if kma_ran else [],
        "status": kma_status
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
