import os

rule nanopore_host_removal:
    input:
        dbjson=DB_JSON,
        fq=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "trimmed_nanopore.fastq.gz")
    output:
        fq=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "nonhost_nanopore.fastq.gz"),
        stats=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "host_stats.json")
    conda:
        os.path.join(ENVS_DIR, "host_nanopore.yaml")
    threads: THREADS
    params:
        preset=lambda wc: config["host_removal"]["minimap2_preset"]
    shell:
        r"""
        python3 scripts/run_minimap2_host_remove.py \
          --dbjson {input.dbjson} \
          --fastq {input.fq} \
          --out_fastq {output.fq} \
          --threads {threads} \
          --preset {params.preset}
        """
