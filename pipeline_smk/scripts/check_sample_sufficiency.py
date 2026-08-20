import argparse
import json
import os
import sys
import gzip

def count_fastq(fastq_path):
    reads = 0
    bases = 0
    open_func = gzip.open if fastq_path.endswith('.gz') else open
    try:
        with open_func(fastq_path, 'rt') as f:
            while True:
                header = f.readline()
                if not header:
                    break
                seq = f.readline().strip()
                f.readline() # +
                f.readline() # qual
                reads += 1
                bases += len(seq)
    except EOFError:
        pass # Handle truncated gzips gracefully if they exist
    return reads, bases

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--r1", help="Path to R1 fastq")
    parser.add_argument("--r2", help="Path to R2 fastq")
    parser.add_argument("--fastq", help="Path to single-end/nanopore fastq")
    parser.add_argument("--min-reads", type=int, default=1000)
    parser.add_argument("--min-bases", type=int, default=100000)
    parser.add_argument("--out", required=True, help="Path to output status JSON")
    
    args = parser.parse_args()

    total_reads = 0
    total_bases = 0
    
    if args.r1 and args.r2:
        r1_reads, r1_bases = count_fastq(args.r1)
        r2_reads, r2_bases = count_fastq(args.r2)
        total_reads = r1_reads + r2_reads
        total_bases = r1_bases + r2_bases
    elif args.fastq:
        total_reads, total_bases = count_fastq(args.fastq)
    else:
        print("Error: Must provide either --r1 and --r2, or --fastq", file=sys.stderr)
        sys.exit(1)
        
    sufficient = (total_reads >= args.min_reads) and (total_bases >= args.min_bases)
    
    status = {
        "sufficient": sufficient,
        "total_reads": total_reads,
        "total_bases": total_bases,
        "min_reads_required": args.min_reads,
        "min_bases_required": args.min_bases
    }
    
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as f:
        json.dump(status, f, indent=4)
        
if __name__ == "__main__":
    main()
