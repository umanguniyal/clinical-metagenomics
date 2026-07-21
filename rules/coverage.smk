import os

rule illumina_coverage_and_mode:
    input:
        contigs=os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "contigs.fasta"),
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R2.fastq.gz")
    output:
        depth=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "contig_depth.txt"),
        mean=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "mean_coverage.txt"),
        mode=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "amr_mode.txt")
    conda:
        os.path.join(ENVS_DIR, "coverage_illumina.yaml")
    threads: THREADS
    shell:
        r"""
        python3 scripts/illumina_coverage_and_mode.py \
          --contigs {input.contigs} --r1 {input.r1} --r2 {input.r2} \
          --out_depth {output.depth} --out_mean {output.mean} --out_mode {output.mode} \
          --threads {threads} --threshold {config[coverage][min_coverage_threshold]}
        """
