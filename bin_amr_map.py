import argparse, csv, json, os

def read_contigs_from_fasta(path):
    contigs = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                contigs.append(line.strip().lstrip(">").split()[0])
    return contigs

def parse_abricate_card(path):
    hits = []
    if not os.path.exists(path):
        return hits
    with open(path) as f:
        header = f.readline()
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 14:
                continue
            # ABRicate: [0]=FILE [1]=SEQUENCE
            contig_id = fields[1].split()[0]
            hits.append({
                "contig": contig_id,
                "gene": fields[5],
                "coverage": fields[8],
                "identity": fields[9],
                "resistance": fields[12]
            })
    return hits

def parse_amrfinder(path):
    hits = []
    if not os.path.exists(path):
        return hits
    with open(path) as f:
        header = f.readline()
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 13:
                continue
            contig_id = fields[1].split()[0]
            hits.append({
                "contig": contig_id,
                "gene": fields[6],
                "subclass": fields[11],
                "element_type": fields[12]
            })
    return hits

def parse_checkm2_quality(path):
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row.get("Name") or row.get("Bin Id") or row.get("bin") or ""
            if not name:
                continue
            out[name] = {
                "completeness": row.get("Completeness", "N/A"),
                "contamination": row.get("Contamination", "N/A"),
            }
    return out

def parse_prokka_organism(annotation_dir):
    org = {}
    if not os.path.isdir(annotation_dir):
        return org
    for bin_name in os.listdir(annotation_dir):
        txt = os.path.join(annotation_dir, bin_name, f"{bin_name}.txt")
        if os.path.exists(txt):
            with open(txt) as f:
                for line in f:
                    if line.startswith("organism:"):
                        org[bin_name] = line.split(":", 1)[1].strip()
                        break
    return org

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin_dir", required=True)
    ap.add_argument("--annotation_dir", required=True)
    ap.add_argument("--contig_abricate_card", required=True)
    ap.add_argument("--contig_amrfinder", required=True)
    ap.add_argument("--checkm2_quality", required=True)
    ap.add_argument("--out_json", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_summary", required=True)
    args = ap.parse_args()

    # contig -> bin
    contig_to_bin = {}
    for fn in sorted(os.listdir(args.bin_dir)):
        if not fn.endswith(".fa"):
            continue
        bin_name = fn.replace(".fa", "")
        for contig in read_contigs_from_fasta(os.path.join(args.bin_dir, fn)):
            contig_to_bin[contig] = bin_name

    bin_quality = parse_checkm2_quality(args.checkm2_quality)
    bin_organism = parse_prokka_organism(args.annotation_dir)

    ab_hits = parse_abricate_card(args.contig_abricate_card)
    amrf_hits = parse_amrfinder(args.contig_amrfinder)

    mapped_ab = []
    for h in ab_hits:
        b = contig_to_bin.get(h["contig"], "unbinned")
        mapped_ab.append({
            **h,
            "bin": b,
            "organism": bin_organism.get(b, "unknown"),
            "bin_completeness": bin_quality.get(b, {}).get("completeness", "N/A"),
            "bin_contamination": bin_quality.get(b, {}).get("contamination", "N/A"),
        })

    mapped_amrf = []
    for h in amrf_hits:
        b = contig_to_bin.get(h["contig"], "unbinned")
        mapped_amrf.append({
            **h,
            "bin": b,
            "organism": bin_organism.get(b, "unknown"),
        })

    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)

    # TSV
    with open(args.out_tsv, "w") as f:
        f.write("Gene\tContig\tBin\tOrganism\tResistance\tIdentity\tCoverage\tCompleteness\tContamination\n")
        for r in mapped_ab:
            f.write(
                f"{r['gene']}\t{r['contig']}\t{r['bin']}\t{r['organism']}\t"
                f"{r.get('resistance','')}\t{r.get('identity','')}\t{r.get('coverage','')}\t"
                f"{r['bin_completeness']}\t{r['bin_contamination']}\n"
            )

    # per-bin summary
    summary = {}
    for r in mapped_ab:
        b = r["bin"]
        summary.setdefault(b, {"organism": r["organism"], "genes": [], "resistances": set()})
        summary[b]["genes"].append(r["gene"])
        if r.get("resistance"):
            summary[b]["resistances"].add(r["resistance"])

    with open(args.out_summary, "w") as f:
        f.write("Bin\tOrganism\tAMR_Gene_Count\tGenes\tResistance_Classes\n")
        for b in sorted(summary.keys()):
            genes = "; ".join(summary[b]["genes"])
            res = "; ".join(sorted(summary[b]["resistances"])) if summary[b]["resistances"] else "N/A"
            f.write(f"{b}\t{summary[b]['organism']}\t{len(summary[b]['genes'])}\t{genes}\t{res}\n")

    report = {
        "total_bins": len(set(contig_to_bin.values())),
        "bins_with_amr": len(set(r["bin"] for r in mapped_ab if r["bin"] != "unbinned")),
        "unbinned_amr_hits": len([r for r in mapped_ab if r["bin"] == "unbinned"]),
        "bin_amr_abricate": mapped_ab,
        "bin_amr_amrfinder": mapped_amrf,
        "bin_quality": bin_quality,
        "bin_organisms": bin_organism
    }

    with open(args.out_json, "w") as f:
        json.dump(report, f, indent=2)

if __name__ == "__main__":
    main()
