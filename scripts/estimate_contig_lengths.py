import argparse

def contig_lengths(fasta):
    lens = []
    l = 0
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if l > 0: lens.append(l)
                l = 0
            else:
                l += len(line.strip())
        if l > 0: lens.append(l)
    return lens

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    lens = contig_lengths(args.fasta)
    if not lens:
        print("0", file=open(args.out, "w"))
        return
    mean = int(sum(lens) / len(lens))
    print(f"{mean}", file=open(args.out, "w"))

if __name__ == "__main__":
    main()
