import os
import argparse
import glob

def parse_benchmark_file(filepath):
    # s, h:m:s, max_rss, max_vms, max_uss, max_pss, io_in, io_out, mean_load, cpu_time
    # Values might be empty or 0.
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
            if len(lines) < 2:
                return None
            headers = lines[0].strip().split('\t')
            values = lines[1].strip().split('\t')
            data = dict(zip(headers, values))
            return data
    except Exception:
        return None

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results_dir", required=True)
    args = parser.parse_args()

    results_dir = os.path.abspath(args.results_dir)
    search_pattern = os.path.join(results_dir, "**", "benchmarks", "*.txt")
    
    # We will use glob.glob with recursive=True
    bench_files = glob.glob(search_pattern, recursive=True)

    out_lines = []
    out_lines.append("Rule\tSample\tPlatform\tWall_Time_Seconds\tPeak_RAM_MB\tIO_Read_MB\tIO_Write_MB")

    total_time = 0.0

    for bf in bench_files:
        # e.g. results_dir / platform / sample / benchmarks / rule.txt
        rel_path = os.path.relpath(bf, results_dir)
        parts = rel_path.split(os.sep)
        
        rule = parts[-1].replace(".txt", "")
        platform = "N/A"
        sample = "N/A"
        
        if len(parts) >= 4 and parts[-2] == "benchmarks":
            sample = parts[-3]
            platform = parts[-4]
            
        data = parse_benchmark_file(bf)
        if data:
            wall_time = float(data.get("s", 0))
            max_rss = float(data.get("max_rss", 0))
            io_in = float(data.get("io_in", 0))
            io_out = float(data.get("io_out", 0))
            
            total_time += wall_time
            
            # max_rss in snakemake is usually in MB, but sometimes it depends on the platform.
            # Assuming max_rss is MB as requested.
            out_lines.append(f"{rule}\t{sample}\t{platform}\t{wall_time:.2f}\t{max_rss:.2f}\t{io_in:.2f}\t{io_out:.2f}")

    out_lines.append(f"TOTAL\t\t\t{total_time:.2f}\t\t\t")

    out_file = "pipeline_benchmark_summary.txt"
    with open(out_file, "w") as f:
        f.write("\n".join(out_lines) + "\n")

if __name__ == "__main__":
    main()
