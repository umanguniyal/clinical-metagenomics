import argparse, csv, json, os, math, re
from collections import Counter

# Bulletproof plotting: gracefully skip if matplotlib/numpy are missing in this conda env
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np
    CAN_PLOT = True
except ImportError:
    CAN_PLOT = False

def calculate_kmer_profile(fasta_path, k=4):
    kmers = Counter()
    total_kmers = 0
    if not os.path.exists(fasta_path): return kmers
    with open(fasta_path) as f:
        seq = "".join([line.strip().upper() for line in f if not line.startswith(">")])
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if "N" not in kmer:
                kmers[kmer] += 1
                total_kmers += 1
    if total_kmers > 0:
        for kmer in kmers: kmers[kmer] /= total_kmers
    return kmers

def calculate_contig_kmer_profiles(fasta_path, k=4):
    profiles = {}
    if not os.path.exists(fasta_path): return profiles
    curr_id, curr_seq = "", []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                if curr_id and curr_seq:
                    seq = "".join(curr_seq).upper()
                    kmers = Counter()
                    tot = 0
                    for i in range(len(seq) - k + 1):
                        kmer = seq[i:i+k]
                        if "N" not in kmer:
                            kmers[kmer] += 1
                            tot += 1
                    if tot > 0:
                        for kmer in kmers: kmers[kmer] /= tot
                    profiles[curr_id] = kmers
                curr_id = line.strip().lstrip(">").split()[0]
                curr_seq = []
            else:
                curr_seq.append(line.strip())
    if curr_id and curr_seq:
        seq = "".join(curr_seq).upper()
        kmers = Counter()
        tot = 0
        for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            if "N" not in kmer:
                kmers[kmer] += 1
                tot += 1
        if tot > 0:
            for kmer in kmers: kmers[kmer] /= tot
        profiles[curr_id] = kmers
    return profiles

def cosine_similarity(prof1, prof2):
    intersection = set(prof1.keys()) & set(prof2.keys())
    numerator = sum([prof1[x] * prof2[x] for x in intersection])
    sum1 = sum([prof1[x]**2 for x in prof1.keys()])
    sum2 = sum([prof2[x]**2 for x in prof2.keys()])
    denominator = math.sqrt(sum1) * math.sqrt(sum2)
    return float(numerator) / denominator if denominator else 0.0

def read_contigs_from_fasta(path):
    contigs = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                contigs.append(line.strip().lstrip(">").split()[0])
    return contigs

def parse_abricate_card(path):
    hits = []
    if not os.path.exists(path): return hits
    with open(path) as f:
        header = f.readline()
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 15: continue
            hits.append({
                "contig": fields[1].split()[0], "gene": fields[5],
                "coverage": fields[9], "identity": fields[10], "resistance": fields[14]
            })
    return hits

def parse_amrfinder(path):
    hits = []
    if not os.path.exists(path): return hits
    with open(path) as f:
        header = f.readline()
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 18: continue
            hits.append({
                "contig": fields[2].split()[0], "gene": fields[6],
                "subclass": fields[12], "element_type": fields[9],
                "coverage": fields[16], "identity": fields[17]
            })
    return hits

def parse_checkm2_quality(path):
    out = {}
    if not os.path.exists(path): return out
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            name = row.get("Name") or row.get("Bin Id") or row.get("bin") or ""
            if name: out[name] = {"completeness": row.get("Completeness", "N/A"), "contamination": row.get("Contamination", "N/A")}
    return out

def parse_centrifuge_bin_taxonomy(tax_dir):
    org = {}
    if not os.path.isdir(tax_dir): return org
    for fn in os.listdir(tax_dir):
        if not fn.endswith(".best_hit.txt"): continue
        bin_name = fn.replace(".best_hit.txt", "")
        with open(os.path.join(tax_dir, fn)) as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 1: 
                    org[bin_name] = parts[0].strip()
    return org

def plot_plasmid_similarities(similarities, outpath):
    if not CAN_PLOT or not similarities:
        return
    plasmids = list(similarities.keys())
    bins = list(next(iter(similarities.values())).keys())
    
    matrix = np.zeros((len(plasmids), len(bins)))
    for i, p in enumerate(plasmids):
        for j, b in enumerate(bins):
            matrix[i, j] = similarities[p].get(b, 0.0)

    fig, ax = plt.subplots(figsize=(max(8, len(bins)), max(4, len(plasmids)*0.5)))
    cax = ax.matshow(matrix, cmap="YlGnBu", vmin=0, vmax=1)
    fig.colorbar(cax, label="K-mer Cosine Similarity")

    ax.set_xticks(range(len(bins)))
    ax.set_xticklabels(bins, rotation=45, ha="left")
    ax.set_yticks(range(len(plasmids)))
    ax.set_yticklabels(plasmids)
    
    for i in range(len(plasmids)):
        for j in range(len(bins)):
            val = matrix[i, j]
            color = "white" if val > 0.6 else "black"
            ax.text(j, i, f"{val:.2f}", ha="center", va="center", color=color, fontsize=8)

    plt.title("Plasmid K-mer Profiler: Similarity to Host Bins", pad=20)
    plt.tight_layout()
    plt.savefig(outpath, dpi=150)
    plt.close()

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin_dir", required=True)
    ap.add_argument("--bin_taxonomy_dir", required=True)
    ap.add_argument("--contig_abricate_card", required=True)
    ap.add_argument("--contig_amrfinder", required=True)
    ap.add_argument("--checkm2_quality", required=True)
    ap.add_argument("--out_json", required=True)
    ap.add_argument("--out_tsv", required=True)
    ap.add_argument("--out_summary", required=True)
    args = ap.parse_args()

    contig_to_bin = {}
    bin_kmer_profiles = {}
    unbinned_fa = None

    for fn in sorted(os.listdir(args.bin_dir)):
        if not fn.endswith(".fa"): continue
        bin_name = fn.replace(".fa", "")
        fa_path = os.path.join(args.bin_dir, fn)
        
        if "unbinned" in bin_name:
            unbinned_fa = fa_path
            for contig in read_contigs_from_fasta(fa_path):
                contig_to_bin[contig] = "unbinned"
        else:
            bin_kmer_profiles[bin_name] = calculate_kmer_profile(fa_path)
            for contig in read_contigs_from_fasta(fa_path):
                contig_to_bin[contig] = bin_name

    similarity_matrix = {}
    if unbinned_fa and bin_kmer_profiles:
        unbinned_profiles = calculate_contig_kmer_profiles(unbinned_fa)
        for contig_id, prof in unbinned_profiles.items():
            if not contig_id.startswith("plasmid"): continue # Only check plasmids
            best_bin, best_sim = None, 0.0
            similarity_matrix[contig_id] = {}
            for bname, bprof in bin_kmer_profiles.items():
                sim = cosine_similarity(prof, bprof)
                similarity_matrix[contig_id][bname] = sim
                if sim > best_sim:
                    best_sim = sim
                    best_bin = bname
            # Rescue unbinned plasmid to bin if DNA accent > 80%
            if best_sim >= 0.80 and best_bin:
                contig_to_bin[contig_id] = best_bin

    # Plot the plasmid graph (will safely do nothing if matplotlib is missing)
    out_graph = os.path.join(os.path.dirname(args.out_json), "plasmid_kmer_similarity.png")
    plot_plasmid_similarities(similarity_matrix, out_graph)

    bin_quality = parse_checkm2_quality(args.checkm2_quality)
    bin_organism = parse_centrifuge_bin_taxonomy(args.bin_taxonomy_dir)
    ab_hits = parse_abricate_card(args.contig_abricate_card)
    amrf_hits = parse_amrfinder(args.contig_amrfinder)

    # Normalize gene names for reliable cross-tool deduplication.
    # ABRicate and AMRFinderPlus can report the same gene with minor differences
    # (trailing whitespace, version suffixes like "_1_AF472622", casing).
    def normalize_gene(name):
        name = name.strip()
        # Strip trailing version/accession suffixes like "_1_AF472622"
        name = re.sub(r'_\d+_[A-Z]{1,2}\d{4,}$', '', name)
        return name.lower()

    dedup_hits = {}
    for h in ab_hits:
        norm_gene = normalize_gene(h["gene"])
        key = (h["contig"], norm_gene)
        b = contig_to_bin.get(h["contig"], "unbinned")
        dedup_hits[key] = {**h, "bin": b, "organism": bin_organism.get(b, "unknown"),
                           "bin_completeness": bin_quality.get(b, {}).get("completeness", "N/A"),
                           "bin_contamination": bin_quality.get(b, {}).get("contamination", "N/A"),
                           "source": ["ABRicate"], "high_confidence": False}

    for h in amrf_hits:
        norm_gene = normalize_gene(h["gene"])
        key = (h["contig"], norm_gene)
        b = contig_to_bin.get(h["contig"], "unbinned")
        if key in dedup_hits:
            if "AMRFinderPlus" not in dedup_hits[key]["source"]: dedup_hits[key]["source"].append("AMRFinderPlus")
            dedup_hits[key]["high_confidence"] = True
        else:
            dedup_hits[key] = {**h, "resistance": h.get("subclass", ""), "bin": b, "organism": bin_organism.get(b, "unknown"),
                               "bin_completeness": bin_quality.get(b, {}).get("completeness", "N/A"),
                               "bin_contamination": bin_quality.get(b, {}).get("contamination", "N/A"),
                               "source": ["AMRFinderPlus"], "high_confidence": False}

    mapped_hits = list(dedup_hits.values())
    os.makedirs(os.path.dirname(args.out_json), exist_ok=True)

    # 1. Output TSV
    with open(args.out_tsv, "w") as f:
        f.write("Gene\tContig\tBin\tOrganism\tResistance\tIdentity\tCoverage\tCompleteness\tContamination\tSources\tHighConfidence\n")
        for r in mapped_hits:
            f.write(f"{r['gene']}\t{r['contig']}\t{r['bin']}\t{r['organism']}\t{r.get('resistance','')}\t{r.get('identity','')}\t{r.get('coverage','')}\t{r['bin_completeness']}\t{r['bin_contamination']}\t{','.join(r['source'])}\t{r['high_confidence']}\n")

    # Metrics for JSON & Summary
    total_bins = len(set(contig_to_bin.values()))
    bins_with_amr = len(set(r["bin"] for r in mapped_hits if r["bin"] != "unbinned"))
    unbinned_amr_hits = len([r for r in mapped_hits if r["bin"] == "unbinned"])

    report = {
        "total_bins": total_bins,
        "bins_with_amr": bins_with_amr,
        "unbinned_amr_hits": unbinned_amr_hits,
        "bin_amr_hits": mapped_hits,
        "bin_quality": bin_quality,
        "bin_organisms": bin_organism,
        "plasmid_graph_path": out_graph if CAN_PLOT else None
    }
    
    # 2. Output JSON
    with open(args.out_json, "w") as f:
        json.dump(report, f, indent=2)
        
    # 3. Output Summary TSV (FIXED: This satisfies Snakemake!)
    with open(args.out_summary, "w") as f:
        f.write("Metric\tValue\n")
        f.write(f"Total_Bins\t{total_bins}\n")
        f.write(f"Bins_With_AMR\t{bins_with_amr}\n")
        f.write(f"Unbinned_AMR_Hits\t{unbinned_amr_hits}\n")

if __name__ == "__main__":
    main()
