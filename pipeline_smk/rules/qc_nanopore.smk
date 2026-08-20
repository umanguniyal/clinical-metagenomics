import os
rule nanopore_nanoplot:
    input:
        fq=lambda wc: NANOPORE_SAMPLES[wc.sample]["fq"]
    output:
        html=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "NanoPlot.html")
    log:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "logs", "nanopore_nanoplot.log")
    benchmark:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "benchmarks", "nanopore_nanoplot.txt")
    conda:
        os.path.join(ENVS_DIR, "nanoqc.yaml")
    threads: THREADS
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p $(dirname {output.html})
        nanoplot -t {threads} --fastq {input.fq} -o $(dirname {output.html})
        """
rule nanopore_nanostat:
    input:
        fq=lambda wc: NANOPORE_SAMPLES[wc.sample]["fq"]
    output:
        stats=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "nanostat.txt")
    log:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "logs", "nanopore_nanostat.log")
    benchmark:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "benchmarks", "nanopore_nanostat.txt")
    conda:
        os.path.join(ENVS_DIR, "nanoqc.yaml")
    threads: 1
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p $(dirname {output.stats})
        nanostat {input.fq} > {output.stats}
        """
rule nanopore_chopper_trim:
    input:
        fq=lambda wc: NANOPORE_SAMPLES[wc.sample]["fq"]
    output:
        fq=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "trimmed_nanopore.fastq.gz")
    log:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "logs", "nanopore_chopper_trim.log")
    benchmark:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "benchmarks", "nanopore_chopper_trim.txt")
    conda:
        os.path.join(ENVS_DIR, "nanoqc.yaml")
    threads: THREADS
    params:
        quality=lambda wc: config["qc"]["nanopore_min_quality"],
        minlength=lambda wc: config["qc"]["nanopore_min_length"],
        maxlength=lambda wc: config["qc"]["nanopore_max_length"],
        trim_approach=lambda wc: config["qc"]["nanopore_trim_approach"],
        headcrop=lambda wc: config["qc"]["nanopore_headcrop"],
        tailcrop=lambda wc: config["qc"]["nanopore_tailcrop"]
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p $(dirname {output.fq})
        python3 scripts/chopper_gz.py \
          --in_fastq {input.fq} \
          --out_fastq {output.fq} \
          --quality {params.quality} \
          --minlength {params.minlength} \
          --maxlength {params.maxlength} \
          --trim_approach {params.trim_approach} \
          --headcrop {params.headcrop} \
          --tailcrop {params.tailcrop} \
          --threads {threads}
        """
