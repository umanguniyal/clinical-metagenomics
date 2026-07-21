#!/usr/bin/env python3
import argparse, os

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)

    rows = []
    print("")
    print("Interactive manifest builder")
    print("Enter samples. Leave platform empty to finish.")
    print("platform must be: illumina or nanopore")
    print("For nanopore: r2 should be blank.")
    print("")

    while True:
        platform = input("platform (illumina/nanopore) [blank to finish]: ").strip()
        if not platform:
            break
        if platform not in ("illumina", "nanopore"):
            print("Invalid platform. Try again.")
            continue

        sample = input("sample_id (e.g. ERR15479252): ").strip()
        if not sample:
            print("sample_id required.")
            continue

        r1 = input("r1 path (.fastq or .fastq.gz): ").strip()
        if not r1:
            print("r1 path required.")
            continue

        r2 = ""
        if platform == "illumina":
            r2 = input("r2 path (.fastq or .fastq.gz): ").strip()
            if not r2:
                print("r2 required for illumina.")
                continue

        rows.append((platform, sample, r1, r2))
        print("Added.\n")

    with open(args.out, "w") as f:
        f.write("platform\tsample\tr1\tr2\n")
        for platform, sample, r1, r2 in rows:
            f.write(f"{platform}\t{sample}\t{r1}\t{r2}\n")

    print(f"\nWrote manifest: {args.out}")
    print(f"Rows: {len(rows)}")

if __name__ == "__main__":
    main()
