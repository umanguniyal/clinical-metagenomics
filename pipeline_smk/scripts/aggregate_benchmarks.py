import sys
import os
import glob
import csv
from datetime import timedelta

def main():
    if len(sys.argv) < 3:
        print("Usage: python aggregate_benchmarks.py <results_dir> <output_txt>")
        sys.exit(1)
        
    results_dir = sys.argv[1]
    out_file = sys.argv[2]
    
    benchmark_files = glob.glob(os.path.join(results_dir, "**", "benchmarks", "*.txt"), recursive=True)
    
    total_time = 0.0
    total_io = 0.0
    
    job_stats = []
    
    for f in benchmark_files:
        try:
            with open(f) as fh:
                reader = csv.DictReader(fh, delimiter='\t')
                for row in reader:
                    # Snakemake benchmark header: s, h:m:s, max_rss, max_vms, max_uss, max_pss, io_in, io_out, mean_load, cpu_time
                    s = float(row.get('s', 0))
                    io_in = float(row.get('io_in', 0))
                    io_out = float(row.get('io_out', 0))
                    
                    job_name = os.path.basename(f).replace('.txt', '')
                    sample_name = os.path.basename(os.path.dirname(os.path.dirname(f)))
                    
                    total_time += s
                    total_io += (io_in + io_out)
                    
                    job_stats.append({
                        'job': job_name,
                        'sample': sample_name,
                        'time_s': s,
                        'io_mb': (io_in + io_out)
                    })
        except Exception as e:
            pass
            
    # Sort by time descending
    job_stats.sort(key=lambda x: x['time_s'], reverse=True)
    
    with open(out_file, 'w') as fh:
        fh.write("=== PIPELINE BENCHMARK REPORT ===\n\n")
        fh.write(f"Total Cumulative Execution Time: {timedelta(seconds=int(total_time))}\n")
        fh.write(f"Total Cumulative I/O (MB): {total_io:.2f}\n\n")
        
        fh.write(f"{'Sample':<20} | {'Job Name':<40} | {'Time (s)':<15} | {'I/O (MB)':<15}\n")
        fh.write("-" * 95 + "\n")
        
        for job in job_stats:
            fh.write(f"{job['sample']:<20} | {job['job']:<40} | {job['time_s']:<15.2f} | {job['io_mb']:<15.2f}\n")

if __name__ == "__main__":
    main()
