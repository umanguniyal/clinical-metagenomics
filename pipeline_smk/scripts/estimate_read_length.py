import argparse, gzip, os, statistics

def open_maybe_gz(path):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", errors="ignore")
    return open(path, "r", errors="ignore")

def sample_lengths_fastq(path, n_reads=2000):
    lens = []
    with open_maybe_gz(path) as f:
        i = 0
        for line_num, line in enumerate(f, start=1):
            if line_num % 4 == 2:  # sequence
                lens.append(len(line.strip()))
                i += 1
                if i >= n_reads:
                    break
    return lens

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--fastq", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--n_reads", type=int, default=2000)
    args = ap.parse_args()

    lens = sample_lengths_fastq(args.fastq, args.n_reads)
    if not lens:
        # fallback
        est = 150
    else:
        est = int(statistics.median(lens))
        if est < 1:
            est = 150

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w") as f:
        f.write(str(est) + "\n")

if __name__ == "__main__":
    main()
