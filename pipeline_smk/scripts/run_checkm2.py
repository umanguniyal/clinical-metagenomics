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

    checkm2_db = db.get("checkm2_db")
    os.makedirs(args.outdir, exist_ok=True)

    # List all files ending in .fa in the bin_dir
    bin_files = sorted(glob.glob(os.path.join(args.bin_dir, "*.fa")))

    # Robustly detect and ignore only "unbinned" bin (the dummy)
    real_bins = [f for f in bin_files if not os.path.basename(f).startswith("bin.unbinned") \
                 and os.path.getsize(f) > 0 and any(line.startswith(">") for line in open(f))]
    
    if not real_bins or not checkm2_db:
        # Write a placeholder quality_report.tsv (just header, i.e. means: nothing to score)
        q = os.path.join(args.outdir, "quality_report.tsv")
        with open(q, "w") as out:
            out.write("Name\tCompleteness\tContamination\n")
        return

    # Let CheckM2 manage the output directory; use --force for reruns.
    cmd = [
        "checkm2", "predict",
        "--input", args.bin_dir,
        "--output-directory", args.outdir,
        "--force",
        "--threads", str(args.threads),
        "--database_path", checkm2_db,
        "-x", "fa",
    ]
    subprocess.check_call(cmd)

if __name__ == "__main__":
    main()
