#!/usr/import/env python3
import argparse
import json
import os
import subprocess
import glob

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dbjson", required=True)
    ap.add_argument("--bin_dir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    with open(args.dbjson) as f:
        db = json.load(f)

    checkv_db = db.get("checkv_db")
    os.makedirs(args.outdir, exist_ok=True)

    bin_files = sorted(glob.glob(os.path.join(args.bin_dir, "*.fa")))
    
    if not bin_files or not checkv_db:
        # Write dummy summary
        q = os.path.join(args.outdir, "quality_summary.tsv")
        with open(q, "w") as out:
            out.write("contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tcompleteness\tcontamination\n")
        return

    # Combine all viral bins into a single fasta for CheckV
    # MetaBAT bins might have overlapping contig names, but usually contig IDs are unique in the assembly.
    combined_fa = os.path.join(args.outdir, "combined_viral_bins.fa")
    with open(combined_fa, "w") as out:
        for b in bin_files:
            with open(b) as f:
                out.write(f.read())
                
    # Run CheckV
    cmd = [
        "checkv", "end_to_end",
        combined_fa,
        args.outdir,
        "-t", str(args.threads),
        "-d", checkv_db
    ]
    subprocess.check_call(cmd)

if __name__ == "__main__":
    main()
