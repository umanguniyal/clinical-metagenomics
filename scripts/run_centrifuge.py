import argparse, json, os, subprocess, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dbjson", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--report", required=True)
    ap.add_argument("--kreport", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--min_hitlen", type=int, default=22)
    ap.add_argument("--r1")
    ap.add_argument("--r2")
    ap.add_argument("--fastq")
    args = ap.parse_args()

    with open(args.dbjson) as f:
        db = json.load(f)
    prefix = db.get("centrifuge_prefix")
    if not prefix:
        # If centrifuge DB missing, create empty outputs (non-fatal)
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        open(args.out, "w").write("SKIPPED\n")
        open(args.report, "w").write("SKIPPED\n")
        open(args.kreport, "w").write("SKIPPED\n")
        return

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    if args.r1 and args.r2:
        cent_input = f"-1 '{args.r1}' -2 '{args.r2}'"
    else:
        cent_input = f"-U '{args.fastq}'"

    cmd = (
        f"centrifuge -x '{prefix}' {cent_input} "
        f"-S '{args.out}' --report-file '{args.report}' "
        f"-p {args.threads} --min-hitlen {args.min_hitlen}"
    )
    subprocess.check_call(["bash", "-lc", cmd])

    kcmd = f"centrifuge-kreport -x '{prefix}' '{args.out}' > '{args.kreport}'"
    subprocess.check_call(["bash", "-lc", kcmd])

if __name__ == "__main__":
    main()
