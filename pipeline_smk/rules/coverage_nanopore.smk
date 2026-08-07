import os

rule nanopore_coverage_and_mode:
    input:
        contigs = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contigs.fasta"),
        fq = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "nonhost_nanopore.fastq.gz")
    output:
        depth = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "contig_depth.txt"),
        mean = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "mean_coverage.txt"),
        mode = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "amr_mode.txt")
    conda:
        os.path.join(ENVS_DIR, "coverage_nanopore.yaml")
    threads: THREADS
    shell:
        r"""
        python3 scripts/nanopore_coverage_and_mode.py \
          --contigs {input.contigs} --fastq {input.fq} \
          --out_depth {output.depth} --out_mean {output.mean} --out_mode {output.mode} \
          --threads {threads} --threshold {config[coverage][min_coverage_threshold]}
        """
