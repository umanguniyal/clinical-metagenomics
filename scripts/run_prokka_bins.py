import argparse, glob, os, subprocess

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bin_dir", required=True)
    ap.add_argument("--out_dir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    bins = sorted(glob.glob(os.path.join(args.bin_dir, "bin.*.fa")))
    if not bins:
        # still ok
        return

    for b in bins:
        name = os.path.basename(b).replace(".fa", "")
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

if __name__ == "__main__":
    main()
