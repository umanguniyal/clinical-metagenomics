import argparse, json, os, subprocess, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dbjson", required=True)
    ap.add_argument("--amr_mode_file", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", default="")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--platform", required=True, choices=["illumina", "nanopore"])
    args = ap.parse_args()

    with open(args.amr_mode_file) as f:
        mode = f.read().strip()

    os.makedirs(args.outdir, exist_ok=True)

    # Skip logic: Skip for both platforms if coverage is high enough for contig AMR.
    # Run read-based KMA only when mode indicates "reads" strategy.
    if mode != "reads_illumina" and mode != "reads_nanopore":
        open(os.path.join(args.outdir, "kma.res"), "w").write("SKIPPED\n")
        open(os.path.join(args.outdir, "kma.tsv"), "w").write("SKIPPED\n")
        return

    with open(args.dbjson) as f:
        db = json.load(f)

    prefix = db.get("resfinder_kma_prefix")
    if not prefix:
        print("ERROR: ResFinder KMA prefix not detected.", file=sys.stderr)
        sys.exit(1)

    outprefix = os.path.join(args.outdir, "kma")

    cmd = ["kma", "-t_db", prefix, "-o", outprefix, "-t", str(args.threads), "-tsv"]
    
    if args.platform == "nanopore":
        cmd += ["-i", args.r1, "-mem_mode"]
    else:
        cmd += ["-ipe", args.r1, args.r2]
    
    subprocess.check_call(cmd)

    res_file = outprefix + ".res"
    tsv_file = outprefix + ".tsv"

    if not os.path.exists(res_file):
        # KMA may not create .res if no hits found; create empty with header
        header = "Template\tIdentity\tQuery_Identity\tQuery_Coverage\tDepth\tExpected_Coverage\tPercent_Expected\n"
        with open(res_file, "w") as f:
            f.write(header)

    final_res = os.path.join(args.outdir, "kma.res")
    
    if os.path.exists(res_file) and res_file != final_res:
        os.rename(res_file, final_res)

    if os.path.exists(tsv_file):
        os.rename(tsv_file, os.path.join(args.outdir, "kma.tsv"))
    else:
        subprocess.check_call(["cp", final_res, os.path.join(args.outdir, "kma.tsv")])

if __name__ == "__main__":
    main()
