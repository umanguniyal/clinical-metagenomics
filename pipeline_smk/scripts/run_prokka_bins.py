import argparse, glob, os, subprocess, json

def get_passing_bins(report, min_completeness=50.0):
    passing_bins = set()
    if not os.path.exists(report):
        return passing_bins
    
    with open(report, "r") as f:
        header = f.readline().strip().split('\t')
        try:
            name_idx = header.index("Name")
            comp_idx = header.index("Completeness")
        except ValueError:
            return passing_bins
            
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > max(name_idx, comp_idx):
                try:
                    if float(parts[comp_idx]) >= min_completeness:
                        passing_bins.add(parts[name_idx])
                except ValueError:
                    pass
    return passing_bins

def parse_prokka_results(out_dir, bin_name):
    txt_file = os.path.join(out_dir, f"{bin_name}.txt")
    tsv_file = os.path.join(out_dir, f"{bin_name}.tsv")
    
    if not os.path.exists(txt_file) or not os.path.exists(tsv_file):
        return None
        
    stats = {
        "cds_count": 0,
        "trna_count": 0,
        "rrna_count": 0,
        "hypothetical_count": 0,
        "hypothetical_ratio": 0.0,
        "gene_density": 0.0,
        "coding_density_pct": 0.0,
        "total_bases": 0
    }
    
    with open(txt_file, "r") as f:
        for line in f:
            if line.startswith("bases:"):
                stats["total_bases"] = int(line.split(":")[1].strip())
            elif line.startswith("CDS:"):
                stats["cds_count"] = int(line.split(":")[1].strip())
            elif line.startswith("tRNA:"):
                stats["trna_count"] = int(line.split(":")[1].strip())
            elif line.startswith("rRNA:"):
                stats["rrna_count"] = int(line.split(":")[1].strip())
                
    total_cds_length = 0
    with open(tsv_file, "r") as f:
        header = f.readline().strip().split('\t')
        try:
            len_idx = header.index("length_bp")
            prod_idx = header.index("product")
            ftype_idx = header.index("ftype")
        except ValueError:
            return stats
            
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) > max(len_idx, prod_idx, ftype_idx):
                if parts[ftype_idx] == "CDS":
                    total_cds_length += int(parts[len_idx])
                    if parts[prod_idx].lower() == "hypothetical protein":
                        stats["hypothetical_count"] += 1
                        
    if stats["cds_count"] > 0:
        stats["hypothetical_ratio"] = round(stats["hypothetical_count"] / stats["cds_count"], 4)
    if stats["total_bases"] > 0:
        stats["gene_density"] = round(stats["cds_count"] / (stats["total_bases"] / 1000000.0), 2)
        stats["coding_density_pct"] = round((total_cds_length / stats["total_bases"]) * 100, 2)
        
    return stats

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin_dir", required=True)
    ap.add_argument("--checkm2_report", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    bins = sorted(glob.glob(os.path.join(args.bin_dir, "bin.*.fa")))
    if not bins:
        with open(os.path.join(args.out_dir, "prokka_stats.json"), "w") as f:
            json.dump({}, f)
        return

    passing_bins = get_passing_bins(args.checkm2_report, 50.0)
    all_stats = {}

    for b in bins:
        name = os.path.basename(b).replace(".fa", "")
        if name in ["bin.tooShort", "bin.unbinned", "bin.lowDepth", "tooShort", "unbinned", "lowDepth"]:
            print(f"Skipping {name}: pseudo-bin")
            continue
        if name not in passing_bins:
            print(f"Skipping {name}: completeness < 50%")
            continue
            
        out = os.path.join(args.out_dir, name)
        os.makedirs(out, exist_ok=True)
        subprocess.call([
            "prokka", b,
            "--outdir", out,
            "--prefix", name,
            "--metagenome",
            "--cpus", str(args.threads),
            "--force",
            "--compliant"
        ])
        
        stats = parse_prokka_results(out, name)
        if stats:
            all_stats[name] = stats
            
    with open(os.path.join(args.out_dir, "prokka_stats.json"), "w") as f:
        json.dump(all_stats, f, indent=2)

if __name__ == "__main__":
    main()
