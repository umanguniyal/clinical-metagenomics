import argparse, os, subprocess

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--contigs", required=True)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--out_depth", required=True)
    ap.add_argument("--out_mean", required=True)
    ap.add_argument("--out_mode", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--threshold", type=float, default=10.0)
    args = ap.parse_args()

    outdir = os.path.dirname(args.out_depth)
    os.makedirs(outdir, exist_ok=True)

    bam = os.path.join(outdir, "contigs_sorted.bam")

    cmd = (
        f"minimap2 -ax map-ont -t {args.threads} '{args.contigs}' '{args.fastq}' | "
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
        mode = "contig_low_nano"

    with open(args.out_mode, "w") as f:
        f.write(mode + "\n")

if __name__ == "__main__":
    main()
