#!/usr/bin/env python3
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

    busco_db = db.get("busco_db")
    os.makedirs(args.outdir, exist_ok=True)

    bin_files = sorted(glob.glob(os.path.join(args.bin_dir, "*.fa")))
    
    if not bin_files or not busco_db:
        # Write dummy summary
        q = os.path.join(args.outdir, "short_summary.txt")
        with open(q, "w") as out:
            out.write("No eukaryotic bins found.\n")
        return

    # Run BUSCO in auto-lineage mode (eukaryota) using the lineages dir
    # Actually, BUSCO downloads lineages automatically if not specified, 
    # but we can set --download_path to busco_db
    for b in bin_files:
        bin_name = os.path.basename(b).replace('.fa', '')
        cmd = [
            "busco", "-i", b,
            "-o", bin_name,
            "--out_path", args.outdir,
            "-m", "genome",
            "-c", str(args.threads),
            "--auto-lineage-euk",
            "--download_path", busco_db,
            "--offline"
        ]
        try:
            subprocess.check_call(cmd)
        except Exception as e:
            print(f"BUSCO failed on {bin_name}: {e}")

    # Combine summaries if any
    summary_file = os.path.join(args.outdir, "short_summary.txt")
    with open(summary_file, "w") as out:
        for b in bin_files:
            bin_name = os.path.basename(b).replace('.fa', '')
            txts = glob.glob(os.path.join(args.outdir, bin_name, f"short_summary.*.txt"))
            if txts:
                with open(txts[0]) as sf:
                    out.write(f"--- {bin_name} ---\n")
                    out.write(sf.read())
                    out.write("\n")

if __name__ == "__main__":
    main()
