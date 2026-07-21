import argparse, gzip, sys, os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--minsize", type=int, default=10000)  # 10 kB minimum
    args = ap.parse_args()
    path = args.fastq
    # 1. Gzip integrity check
    try:
        if path.endswith(".gz"):
            f = gzip.open(path, "rt", errors="ignore")
        else:
            f = open(path, "rt", errors="ignore")
    except Exception as e:
        print(f"FAIL: Cannot open/read (gz error): {e}", file=sys.stderr)
        sys.exit(2)
    # 2. Check for at least one record starting with @
    found = False
    n = 0
    for line in f:
        if line.startswith('@'):
            found = True
            break
        n += 1
        if n > 5:
            break
    f.close()
    if not found:
        print("FAIL: No record starts with '@'", file=sys.stderr)
        sys.exit(3)
    # 3. Check min size
    if os.path.getsize(path) < args.minsize:
        print(f"FAIL: File is very small (<{args.minsize} bytes)", file=sys.stderr)
        sys.exit(4)
    sys.exit(0)

if __name__ == "__main__":
    main()
