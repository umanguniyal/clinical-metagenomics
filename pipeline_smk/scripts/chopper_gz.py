import argparse, gzip, subprocess, sys, os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in_fastq", required=True)
    ap.add_argument("--out_fastq", required=True)
    ap.add_argument("--quality", type=int, required=True)
    ap.add_argument("--minlength", type=int, required=True)
    ap.add_argument("--maxlength", type=int, required=True)
    ap.add_argument("--headcrop", type=int, required=True)
    ap.add_argument("--tailcrop", type=int, required=True)
    ap.add_argument("--threads", type=int, default=8)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_fastq), exist_ok=True)

    # Stream input (gz or plain) -> chopper -> gzip output
    if args.in_fastq.endswith(".gz"):
        in_cmd = ["bash", "-lc", f"zcat '{args.in_fastq}'"]
    else:
        in_cmd = ["bash", "-lc", f"cat '{args.in_fastq}'"]

    chopper_cmd = [
        "chopper",
        "--quality", str(args.quality),
        "--minlength", str(args.minlength),
        "--maxlength", str(args.maxlength),
        "--headcrop", str(args.headcrop),
        "--tailcrop", str(args.tailcrop),
        "--threads", str(args.threads),
    ]

    p1 = subprocess.Popen(in_cmd, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(chopper_cmd, stdin=p1.stdout, stdout=subprocess.PIPE)
    p1.stdout.close()

    with gzip.open(args.out_fastq, "wb") as out_gz:
        for chunk in iter(lambda: p2.stdout.read(1024 * 1024), b""):
            out_gz.write(chunk)

    rc = p2.wait()
    if rc != 0:
        sys.exit(rc)

if __name__ == "__main__":
    main()
