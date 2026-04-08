import argparse, json, os, subprocess, sys, glob

def find_available_kmer_sizes(kraken_db):
    """Find all available k-mer distribution files."""
    pattern = os.path.join(kraken_db, "database*mers.kmer_distrib")
    files = glob.glob(pattern)
    
    sizes = []
    for f in files:
        basename = os.path.basename(f)
        try:
            size = int(basename.replace("database", "").replace("mers.kmer_distrib", ""))
            sizes.append(size)
        except:
            pass
    
    return sorted(sizes, reverse=True)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dbjson", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("--readlen_file", required=True)
    ap.add_argument("--kraken_report", required=True)
    ap.add_argument("--out_report", required=True)
    ap.add_argument("--out_species", required=True)
    ap.add_argument("--level", default="S")
    args = ap.parse_args()

    with open(args.dbjson) as f:
        db = json.load(f)

    kraken_db = db.get("kraken_db")
    if not kraken_db:
        print("ERROR: Kraken DB not detected.", file=sys.stderr)
        sys.exit(1)

    try:
        with open(args.readlen_file) as f:
            readlen = int(f.read().strip())
    except:
        readlen = 150 if args.platform == "illumina" else 1000

    outdir = os.path.dirname(args.out_report)
    os.makedirs(outdir, exist_ok=True)

    # Check available k-mers
    available_kmers = find_available_kmer_sizes(kraken_db)
    if not available_kmers:
        print(f"ERROR: No k-mer distribution files found in {kraken_db}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Platform: {args.platform}, Readlen: {readlen}bp", file=sys.stderr)
    print(f"Available k-mers: {available_kmers}", file=sys.stderr)

    # Bracken outputs to the base name (no extension)
    bracken_output_base = args.out_report.replace(".tsv", "")
    
    cmd = [
        "bracken",
        "-d", kraken_db,
        "-i", args.kraken_report,
        "-o", bracken_output_base,
        "-r", str(readlen),
        "-l", args.level
    ]
    
    print(f"Running: {' '.join(cmd)}", file=sys.stderr)
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    print(f"Bracken return code: {result.returncode}", file=sys.stderr)
    
    if result.returncode != 0:
        print(f"Bracken stderr: {result.stderr}", file=sys.stderr)
        print(f"Bracken failed - using Kraken2 report as fallback...", file=sys.stderr)
        
        if os.path.exists(args.kraken_report):
            subprocess.check_call(["cp", args.kraken_report, args.out_report])
            subprocess.check_call(["cp", args.kraken_report, args.out_species])
            print(f"Using Kraken2 classifications as fallback", file=sys.stderr)
        else:
            print(f"ERROR: Kraken2 report not found: {args.kraken_report}", file=sys.stderr)
            sys.exit(1)
    else:
        print(f"Bracken completed successfully", file=sys.stderr)
        
        # Bracken creates output without extension, rename to .tsv
        if os.path.exists(bracken_output_base):
            subprocess.check_call(["mv", bracken_output_base, args.out_report])
            subprocess.check_call(["cp", args.out_report, args.out_species])
            print(f"Moved Bracken output to: {args.out_report}", file=sys.stderr)
        else:
            print(f"ERROR: Bracken output not found at {bracken_output_base}", file=sys.stderr)
            sys.exit(1)

if __name__ == "__main__":
    main()
