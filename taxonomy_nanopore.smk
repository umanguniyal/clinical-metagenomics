import os

rule nanopore_kraken2:
    input:
        dbjson=DB_JSON,
        fq=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "nonhost_nanopore.fastq.gz")
    output:
        kraken_out=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "kraken2_output.txt"),
        kraken_report=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "kraken2_report.txt")
    conda:
        os.path.join(ENVS_DIR, "kraken.yaml")
    threads: THREADS
    params:
        confidence=lambda wc: config["taxonomy"]["kraken2_confidence"]
    shell:
        r"""
        python3 scripts/run_kraken2.py \
          --dbjson {input.dbjson} \
          --out {output.kraken_out} \
          --report {output.kraken_report} \
          --threads {threads} \
          --confidence {params.confidence} \
          --use_names \
          --fastq {input.fq}
        """

rule nanopore_tb_flag:
    input:
        kraken_report=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "kraken2_report.txt")
    output:
        flag=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "tb_flag.txt")
    params:
        threshold=lambda wc: config["taxonomy"]["tb_threshold_percent"]
    shell:
        r"""
        python3 scripts/tb_flag_from_kraken_report.py \
          --report {input.kraken_report} \
          --out {output.flag} \
          --threshold {params.threshold}
        """

rule nanopore_tbprofiler:
    input:
        flag=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "tb_flag.txt"),
        fq=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "nonhost_nanopore.fastq.gz")
    output:
        results=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "tbprofiler", "results", "{sample}.results.json")
    conda:
        os.path.join(ENVS_DIR, "tbprofiler.yaml")
    threads: THREADS
    shell:
        r"""
        python3 scripts/run_tbprofiler.py \
          --flag {input.flag} \
          --platform nanopore \
          --sample {wildcards.sample} \
          --outdir $(dirname $(dirname {output.results})) \
          --threads {threads} \
          --fastq {input.fq}
        """

rule nanopore_centrifuge:
    input:
        dbjson=DB_JSON,
        fq=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "nonhost_nanopore.fastq.gz")
    output:
        centrifuge_out=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "centrifuge_output.txt"),
        centrifuge_report=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "centrifuge_report.tsv"),
        centrifuge_kreport=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "centrifuge_kreport.txt")
    conda:
        os.path.join(ENVS_DIR, "centrifuge.yaml")
    threads: THREADS
    params:
        min_hitlen=lambda wc: config["taxonomy"]["centrifuge_min_hitlen_nanopore"]
    shell:
        r"""
        python3 scripts/run_centrifuge.py \
          --dbjson {input.dbjson} \
          --out {output.centrifuge_out} \
          --report {output.centrifuge_report} \
          --kreport {output.centrifuge_kreport} \
          --threads {threads} \
          --min_hitlen {params.min_hitlen} \
          --fastq {input.fq}
        """

rule nanopore_krona_kraken2:
    input:
        kraken_report=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "kraken2_report.txt"),
        dbjson=DB_JSON
    output:
        html=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "krona_kraken2.html")
    conda:
        os.path.join(ENVS_DIR, "krona.yaml")
    shell:
        r"""
        mkdir -p $(dirname {output.html})
        ktImportTaxonomy -m 3 -o {output.html} {input.kraken_report}
        """

rule nanopore_krona_centrifuge:
    input:
        centrifuge_kreport=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "centrifuge_kreport.txt"),
        dbjson=DB_JSON
    output:
        html=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "krona_centrifuge.html")
    conda:
        os.path.join(ENVS_DIR, "krona.yaml")
    shell:
        r"""
        mkdir -p $(dirname {output.html})
        ktImportTaxonomy -m 3 -o {output.html} {input.centrifuge_kreport}
        """
