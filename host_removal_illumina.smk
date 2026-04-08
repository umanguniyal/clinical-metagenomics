import os

rule illumina_host_removal:
    input:
        dbjson=DB_JSON,
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "trimmed_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "trimmed_R2.fastq.gz")
    output:
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R2.fastq.gz")
    conda:
        os.path.join(ENVS_DIR, "bowtie2.yaml")
    threads: THREADS
    params:
        sensitivity=lambda wc: config["host_removal"]["bowtie2_sensitivity"]
    shell:
        r"""
        python3 scripts/run_bowtie2_host_remove.py \
          --dbjson {input.dbjson} \
          --r1 {input.r1} \
          --r2 {input.r2} \
          --out_r1 {output.r1} \
          --out_r2 {output.r2} \
          --threads {threads} \
          --sensitivity {params.sensitivity}
        """
