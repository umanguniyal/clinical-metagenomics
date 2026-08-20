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

    # Let's run CheckM2 individually on each real bin so that one DIAMOND failure
    # (e.g. on a non-bacterial bin) doesn't crash the entire run.
    tmp_out = os.path.join(args.outdir, "tmp")
    os.makedirs(tmp_out, exist_ok=True)
    
    successful_bins = []
    
    for bin_file in real_bins:
        bin_name = os.path.basename(bin_file).replace(".fa", "")
        bin_tmp = os.path.join(tmp_out, bin_name)
        
        cmd = [
            "checkm2", "predict",
            "--input", bin_file,
            "--output-directory", bin_tmp,
            "--force",
            "--threads", str(args.threads),
            "--database_path", checkm2_db,
            "-x", "fa",
        ]
        
        try:
            subprocess.check_call(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            report_file = os.path.join(bin_tmp, "quality_report.tsv")
            if os.path.exists(report_file):
                successful_bins.append((bin_name, report_file))
        except subprocess.CalledProcessError:
            print(f"CheckM2 failed on {bin_name}, skipping.")
            continue
            
    # Merge reports
    final_report = os.path.join(args.outdir, "quality_report.tsv")
    with open(final_report, "w") as out:
        out.write("Name\tCompleteness\tContamination\n")
        for bin_name, report_file in successful_bins:
            with open(report_file, "r") as f:
                lines = f.readlines()
                if len(lines) > 1:
                    out.write(lines[1]) # The data line
                    
    # Clean up tmp
    import shutil
    shutil.rmtree(tmp_out)

if __name__ == "__main__":
    main()
