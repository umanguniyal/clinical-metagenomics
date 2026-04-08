import argparse, gzip, json, os, subprocess, sys, tempfile

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dbjson", required=True)
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--out_fastq", required=True)
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--preset", default="map-ont")
    args = ap.parse_args()

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

if __name__ == "__main__":
    main()
