import os
import glob
import argparse

def check_milestone(sample_dir, file_rel_path):
    path = os.path.join(sample_dir, file_rel_path)
    if os.path.exists(path) and os.path.getsize(path) > 0:
        return "✓"
    return "✗"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results_dir", required=True)
    args = parser.parse_args()

    results_dir = args.results_dir
    print(f"\nScanning results directory: {results_dir}\n")
    
    print(f"{'Sample':<20} {'Platform':<10} {'QC':<5} {'Assembly':<10} {'Coverage':<10} {'AMR':<5} {'Binning':<10} {'Report':<8} {'Status'}")
    print("-" * 90)

    for platform in ["illumina", "nanopore"]:
        platform_dir = os.path.join(results_dir, platform)
        if not os.path.isdir(platform_dir):
            continue
        
        for sample in os.listdir(platform_dir):
            sample_dir = os.path.join(platform_dir, sample)
            if not os.path.isdir(sample_dir):
                continue
                
            qc = check_milestone(sample_dir, "qc/fastp.json") if platform == "illumina" else check_milestone(sample_dir, "qc/trimmed_nanopore.fastq.gz")
            assembly = check_milestone(sample_dir, "assembly/assembly.fasta")
            coverage = check_milestone(sample_dir, "coverage/coverage_mode.json")
            amr = check_milestone(sample_dir, "amr/abricate_card.tsv")
            binning = check_milestone(sample_dir, "bin_taxonomy/bin_organisms.json")
            report = check_milestone(sample_dir, "report/report_bioinformatician.json")
            
            # Simple status determination
            status = "COMPLETE"
            if report == "✗":
                status = "FAILED"
                if binning == "✗": status = "FAILED at: binning/taxonomy"
                if amr == "✗": status = "FAILED at: amr"
                if coverage == "✗": status = "FAILED at: coverage"
                if assembly == "✗": status = "FAILED at: assembly"
                if qc == "✗": status = "FAILED at: qc"
                
            print(f"{sample:<20} {platform:<10} {qc:<5} {assembly:<10} {coverage:<10} {amr:<5} {binning:<10} {report:<8} {status}")

if __name__ == "__main__":
    main()
