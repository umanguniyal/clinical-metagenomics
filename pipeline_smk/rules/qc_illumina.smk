import os



rule illumina_fastp_trim:
    input:
        r1=lambda wc: ILLUMINA_SAMPLES[wc.sample]["r1"],
        r2=lambda wc: ILLUMINA_SAMPLES[wc.sample]["r2"]
    output:
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "trimmed_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "trimmed_R2.fastq.gz"),
        json=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "fastp.json"),
        html=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "fastp.html")
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_fastp_trim.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_fastp_trim.txt")
    conda:
        os.path.join(ENVS_DIR, "fastp.yaml")
    threads: THREADS
    params:
        quality=lambda wc: config["qc"]["illumina_min_quality"],
        min_length=lambda wc: config["qc"]["illumina_min_length"],
        window_size=lambda wc: config["qc"]["illumina_cut_window_size"],
        mean_quality=lambda wc: config["qc"]["illumina_cut_mean_quality"]
    shell:
        r"""
        mkdir -p $(dirname {output.r1})
        fastp \
          -i {input.r1} -I {input.r2} \
          -o {output.r1} -O {output.r2} \
          -q {params.quality} \
          -l {params.min_length} \
          -W {params.window_size} \
          -M {params.mean_quality} \
          --json {output.json} \
          --html {output.html} \
          --thread {threads}
        """
