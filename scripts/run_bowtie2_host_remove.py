import argparse, json, os, subprocess, sys

def count_reads(fq_path):
    cmd = f"zcat '{fq_path}' | wc -l"
    try:
        lines = int(subprocess.check_output(["bash", "-c", cmd]).strip())
        return lines // 4
    except:
        return 0

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dbjson", required=True)
    ap.add_argument("--r1", required=True)
    ap.add_argument("--r2", required=True)
    ap.add_argument("--out_r1", required=True)
    ap.add_argument("--out_r2", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--sensitivity", default="very-sensitive")
    args = ap.parse_args()

    total_r1_reads = count_reads(args.r1)

    with open(args.dbjson) as f:
        db = json.load(f)
    host_prefix = db.get("host_bt2_prefix")
    if not host_prefix:
        print("ERROR: Bowtie2 host index prefix not detected (*.1.bt2).", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.out_r1), exist_ok=True)

    # bowtie2 writes nonhost as nonhost_R1/nonhost_R2 using --un-conc
    out_prefix = os.path.join(os.path.dirname(args.out_r1), "nonhost_R%.fastq.gz")

    cmd = [
        "bowtie2",
        "-p", str(args.threads),
        "-x", host_prefix,
        "-1", args.r1,
        "-2", args.r2,
        f"--{args.sensitivity}",
        "--un-conc-gz", out_prefix,
        "-S", os.path.join(os.path.dirname(args.out_r1), "host_mapped.sam"),
    ]
    subprocess.check_call(cmd)

    # ensure expected names
    # bowtie2 outputs nonhost_R1.fastq.gz and nonhost_R2.fastq.gz
    produced_r1 = os.path.join(os.path.dirname(args.out_r1), "nonhost_R1.fastq.gz")
    produced_r2 = os.path.join(os.path.dirname(args.out_r2), "nonhost_R2.fastq.gz")

    if not os.path.exists(produced_r1) or not os.path.exists(produced_r2):
        print("ERROR: bowtie2 did not produce expected nonhost files.", file=sys.stderr)
        sys.exit(1)

    os.rename(produced_r1, args.out_r1)
    os.rename(produced_r2, args.out_r2)

    sam = os.path.join(os.path.dirname(args.out_r1), "host_mapped.sam")
    if os.path.exists(sam):
        os.remove(sam)

    nonhost_r1_reads = count_reads(args.out_r1)
    total_reads = total_r1_reads * 2
    nonhost_reads = nonhost_r1_reads * 2
    host_reads = total_reads - nonhost_reads
    host_percent = (host_reads / total_reads * 100) if total_reads > 0 else 0

    stats = {
        "total_reads_before": total_reads,
        "total_reads_after": nonhost_reads,
        "human_reads_removed": host_reads,
        "human_percent": round(host_percent, 2)
    }

    stats_file = os.path.join(os.path.dirname(args.out_r1), "host_stats.json")
    with open(stats_file, "w") as f:
        json.dump(stats, f, indent=2)

if __name__ == "__main__":
    main()
