import argparse, os, subprocess, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--out_depth", required=True)
    ap.add_argument("--out_mean", required=True)
    ap.add_argument("--out_mode", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--threshold", type=float, default=10.0)
    args = ap.parse_args()

    outdir = os.path.dirname(args.out_depth)
    os.makedirs(outdir, exist_ok=True)

    idx_prefix = os.path.join(outdir, "contigs_index")
    bam = os.path.join(outdir, "contigs_sorted.bam")

    subprocess.check_call(["bowtie2-build", "--threads", str(args.threads), args.contigs, idx_prefix])

    cmd = (
        f"bowtie2 -p {args.threads} -x '{idx_prefix}' -1 '{args.r1}' -2 '{args.r2}' "
        f"2> '{outdir}/bowtie2_map.log' | "
        f"samtools view -@ {args.threads} -bS - | "
        f"samtools sort -@ {args.threads} -o '{bam}'"
    )
    subprocess.check_call(["bash", "-lc", cmd])
    subprocess.check_call(["samtools", "index", bam])

    subprocess.check_call([
        "jgi_summarize_bam_contig_depths",
        "--outputDepth", args.out_depth,
        bam
    ])

    # compute mean coverage from depth.txt col3 (same as your bash)
    mean = 0.0
    n = 0
    with open(args.out_depth) as f:
        next(f, None)
        for line in f:
            parts = line.strip().split()
            if len(parts) >= 3:
                try:
                    mean += float(parts[2])
                    n += 1
                except:
                    pass
    mean_cov = (mean / n) if n else 0.0

    with open(args.out_mean, "w") as f:
        f.write(f"{mean_cov:.2f}\n")

    if mean_cov >= args.threshold:
        mode = "contig_high"
    else:
        mode = "reads_illumina"

    with open(args.out_mode, "w") as f:
        f.write(mode + "\n")

if __name__ == "__main__":
    main()
