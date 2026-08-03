import argparse, json, subprocess, sys, os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dbjson", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--confidence", type=float, default=0.0)
    ap.add_argument("--use_names", action="store_true")
    ap.add_argument("--paired", action="store_true")
    ap.add_argument("--r1")
    ap.add_argument("--r2")
    ap.add_argument("--fastq")
    args = ap.parse_args()

    with open(args.dbjson) as f:
        db = json.load(f)

    kraken_db = db.get("kraken_db")
    if not kraken_db:
        print("ERROR: Kraken DB not detected.", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    cmd = ["kraken2", "--db", kraken_db,
           "--output", args.out,
           "--report", args.report,
           "--threads", str(args.threads)]

    if args.use_names:
        cmd.append("--use-names")

    if args.confidence and args.confidence > 0:
        cmd += ["--confidence", str(args.confidence)]

    if args.paired:
        if not args.r1 or not args.r2:
            print("ERROR: paired mode requires --r1 and --r2", file=sys.stderr)
            sys.exit(1)
        cmd += ["--paired", args.r1, args.r2]
    else:
        if not args.fastq:
            print("ERROR: single mode requires --fastq", file=sys.stderr)
            sys.exit(1)
        cmd.append(args.fastq)

    subprocess.check_call(cmd)

if __name__ == "__main__":
    main()
