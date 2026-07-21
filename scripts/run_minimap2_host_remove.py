import argparse, gzip, json, os, subprocess, sys, tempfile

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
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--out_fastq", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--preset", default="map-ont")
    args = ap.parse_args()

    total_reads = count_reads(args.fastq)

    with open(args.dbjson) as f:
        db = json.load(f)

    host_mmi = db.get("host_mmi")
    if not host_mmi:
        print("ERROR: host MMI not detected (GRCh38*.mmi).", file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.out_fastq), exist_ok=True)

    # minimap2 -> samtools view -f 4 -> samtools fastq -> gzip
    # Need samtools in env as well (host_nanopore.yaml includes it)
    cmd = (
        f"minimap2 -ax {args.preset} -t {args.threads} --secondary=no '{host_mmi}' '{args.fastq}' "
        f"| samtools view -@ {args.threads} -bS -f 4 - "
        f"| samtools fastq -@ {args.threads} - "
        f"| gzip -c > '{args.out_fastq}'"
    )
    subprocess.check_call(["bash", "-lc", cmd])

    nonhost_reads = count_reads(args.out_fastq)
    host_reads = total_reads - nonhost_reads
    host_percent = (host_reads / total_reads * 100) if total_reads > 0 else 0

    stats = {
        "total_reads_before": total_reads,
        "total_reads_after": nonhost_reads,
        "human_reads_removed": host_reads,
        "human_percent": round(host_percent, 2)
    }

    stats_file = os.path.join(os.path.dirname(args.out_fastq), "host_stats.json")
    with open(stats_file, "w") as f:
        json.dump(stats, f, indent=2)

if __name__ == "__main__":
    main()
