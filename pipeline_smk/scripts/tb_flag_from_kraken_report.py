import argparse, os, sys

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--threshold", type=float, default=1.0)
    args = ap.parse_args()

    tb_pct = 0.0
    if os.path.exists(args.report):
        with open(args.report) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 6:
                    pct, rank, name = parts[0], parts[3].strip(), parts[5].strip()
                    if rank == "S" and "Mycobacterium tuberculosis" in name:
                        try:
                            tb_pct = float(pct)
                        except:
                            tb_pct = 0.0

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    if tb_pct > args.threshold:
        open(args.out, "w").write(f"TB_DETECTED\t{tb_pct}\n")
    else:
        # still write file for snakemake; empty means "not detected"
        open(args.out, "w").write("")

if __name__ == "__main__":
    main()
