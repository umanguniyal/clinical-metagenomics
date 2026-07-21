import argparse, os, subprocess

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--minid", type=int, default=80)
    ap.add_argument("--mincov", type=int, default=60)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # ABRicate: card ncbi resfinder vfdb
    for db in ["card", "ncbi", "resfinder", "vfdb"]:
        out = os.path.join(args.outdir, f"abricate_{db}.tsv")
        cmd = [
            "abricate", "--db", db,
            "--threads", str(args.threads),
            "--minid", str(args.minid),
            "--mincov", str(args.mincov),
            args.contigs
        ]
        with open(out, "w") as f:
            subprocess.check_call(cmd, stdout=f)

    # summary
    summ = os.path.join(args.outdir, "abricate_summary.tsv")
    subprocess.check_call(["bash", "-lc", f"abricate --summary {args.outdir}/abricate_*.tsv > {summ}"])

    # AMRFinderPlus
    amrf = os.path.join(args.outdir, "amrfinder_results.tsv")
    subprocess.check_call([
        "amrfinder",
        "-n", args.contigs,
        "-o", amrf,
        "--threads", str(args.threads),
        "--plus",
        "--name", "sample"
    ])

if __name__ == "__main__":
    main()
