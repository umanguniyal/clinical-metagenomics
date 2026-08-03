import argparse, os, subprocess

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--flag", required=True)
    ap.add_argument("--platform", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--r1")
    ap.add_argument("--r2")
    ap.add_argument("--fastq")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # if flag file empty -> skip
    if os.path.exists(args.flag) and os.path.getsize(args.flag) == 0:
        # create results placeholder
        resdir = os.path.join(args.outdir, "results")
        os.makedirs(resdir, exist_ok=True)
        open(os.path.join(resdir, f"{args.sample}.results.json"), "w").write("{}\n")
        return

    if args.platform == "illumina":
        cmd = [
            "tb-profiler", "profile",
            "-1", args.r1,
            "-2", args.r2,
            "-t", str(args.threads),
            "-p", args.sample,
            "--dir", args.outdir,
            "--txt"
        ]
    else:
        cmd = [
            "tb-profiler", "profile",
            "-1", args.fastq,
            "-t", str(args.threads),
            "-p", args.sample,
            "--dir", args.outdir,
            "--platform", "nanopore",
            "--txt"
        ]

    subprocess.check_call(cmd)

if __name__ == "__main__":
    main()
