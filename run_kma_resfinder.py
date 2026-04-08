import argparse, json, os, subprocess, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dbjson", required=True)
    ap.add_argument("--amr_mode_file", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    with open(args.amr_mode_file) as f:
        mode = f.read().strip()

    os.makedirs(args.outdir, exist_ok=True)

    if mode != "reads_illumina":
        open(os.path.join(args.outdir, "kma.res"), "w").write("SKIPPED\n")
        open(os.path.join(args.outdir, "kma.tsv"), "w").write("SKIPPED\n")
        return

    with open(args.dbjson) as f:
        db = json.load(f)

    prefix = db.get("resfinder_kma_prefix")
    if not prefix:
        print("ERROR: ResFinder KMA prefix not detected (resfinder_db/resfinder_db_kma.seq.b).", file=sys.stderr)
        sys.exit(1)

    outprefix = os.path.join(args.outdir, "kma")

    cmd = [
        "kma",
        "-ipe", args.r1, args.r2,
        "-t_db", prefix,
        "-o", outprefix,
        "-t", str(args.threads),
        "-tsv"
    ]
    subprocess.check_call(cmd)

    res_file = outprefix + ".res"
    tsv_file = outprefix + ".tsv"

    if not os.path.exists(res_file):
        print("ERROR: KMA did not create .res output", file=sys.stderr)
        sys.exit(1)

    os.rename(res_file, os.path.join(args.outdir, "kma.res"))

    if os.path.exists(tsv_file):
        os.rename(tsv_file, os.path.join(args.outdir, "kma.tsv"))
    else:
        # fallback: still produce something
        subprocess.check_call(["cp", os.path.join(args.outdir, "kma.res"), os.path.join(args.outdir, "kma.tsv")])

if __name__ == "__main__":
    main()
