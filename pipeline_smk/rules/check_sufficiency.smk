import os
import json
rule check_illumina_sufficiency:
    input:
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R2.fastq.gz")
    output:
        status=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "sample_status.json")
    params:
        min_reads=lambda wc: config.get("sufficiency", {}).get("min_reads_after_host_removal", 1000),
        min_bases=lambda wc: config.get("sufficiency", {}).get("min_bases_after_host_removal", 100000)
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "check_sufficiency.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "check_illumina_sufficiency.txt")
    conda:
        os.path.join(ENVS_DIR, "base_utils.yaml")
    shell:
        r"""
        exec > {log} 2>&1
        python3 scripts/check_sample_sufficiency.py \
            --r1 {input.r1} \
            --r2 {input.r2} \
            --min-reads {params.min_reads} \
            --min-bases {params.min_bases} \
            --out {output.status} > {log} 2>&1
        """
rule check_nanopore_sufficiency:
    input:
        fq=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "nonhost_nanopore.fastq.gz")
    output:
        status=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "sample_status.json")
    params:
        min_reads=lambda wc: config.get("sufficiency", {}).get("min_reads_after_host_removal", 500), # Nanopore usually has fewer but longer reads
        min_bases=lambda wc: config.get("sufficiency", {}).get("min_bases_after_host_removal", 100000)
    log:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "logs", "check_sufficiency.log")
    benchmark:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "benchmarks", "check_nanopore_sufficiency.txt")
    conda:
        os.path.join(ENVS_DIR, "base_utils.yaml")
    shell:
        r"""
        exec > {log} 2>&1
        python3 scripts/check_sample_sufficiency.py \
            --fastq {input.fq} \
            --min-reads {params.min_reads} \
            --min-bases {params.min_bases} \
            --out {output.status} > {log} 2>&1
        """
