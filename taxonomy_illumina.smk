import os

rule illumina_kraken2:
    input:
        dbjson=DB_JSON,
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R2.fastq.gz")
    output:
        kraken_out=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "kraken2_output.txt"),
        kraken_report=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "kraken2_report.txt")
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
          --paired \
          --r1 {input.r1} --r2 {input.r2}
        """

rule illumina_tb_flag:
    input:
        kraken_report=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "kraken2_report.txt")
    output:
        flag=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "tb_flag.txt")
    params:
        threshold=lambda wc: config["taxonomy"]["tb_threshold_percent"]
    shell:
        r"""
        python3 scripts/tb_flag_from_kraken_report.py \
          --report {input.kraken_report} \
          --out {output.flag} \
          --threshold {params.threshold}
        """

rule illumina_tbprofiler:
    input:
        flag=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "tb_flag.txt"),
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R2.fastq.gz")
    output:
        results=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "tbprofiler", "results", "{sample}.results.json")
    conda:
        os.path.join(ENVS_DIR, "tbprofiler.yaml")
    threads: THREADS
    shell:
        r"""
        python3 scripts/run_tbprofiler.py \
          --flag {input.flag} \
          --platform illumina \
          --sample {wildcards.sample} \
          --outdir $(dirname $(dirname {output.results})) \
          --threads {threads} \
          --r1 {input.r1} --r2 {input.r2}
        """

rule illumina_centrifuge:
    input:
        dbjson=DB_JSON,
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R2.fastq.gz")
    output:
        centrifuge_out=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_output.txt"),
        centrifuge_report=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_report.tsv"),
        centrifuge_kreport=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_kreport.txt")
    conda:
        os.path.join(ENVS_DIR, "centrifuge.yaml")
    threads: THREADS
    params:
        min_hitlen=lambda wc: config["taxonomy"]["centrifuge_min_hitlen_illumina"]
    shell:
        r"""
        python3 scripts/run_centrifuge.py \
          --dbjson {input.dbjson} \
          --out {output.centrifuge_out} \
          --report {output.centrifuge_report} \
          --kreport {output.centrifuge_kreport} \
          --threads {threads} \
          --min_hitlen {params.min_hitlen} \
          --r1 {input.r1} --r2 {input.r2}
        """

rule illumina_krona_kraken2:
    input:
        kraken_report=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "kraken2_report.txt"),
        dbjson=DB_JSON
    output:
        html=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "krona_kraken2.html")
    conda:
        os.path.join(ENVS_DIR, "krona.yaml")
    shell:
        r"""
        mkdir -p $(dirname {output.html})
        ktImportTaxonomy -m 3 -o {output.html} {input.kraken_report}
        """

rule illumina_krona_centrifuge:
    input:
        centrifuge_kreport=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_kreport.txt"),
        dbjson=DB_JSON
    output:
        html=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "krona_centrifuge.html")
    conda:
        os.path.join(ENVS_DIR, "krona.yaml")
    shell:
        r"""
        mkdir -p $(dirname {output.html})
        ktImportTaxonomy -m 3 -o {output.html} {input.centrifuge_kreport}
        """
