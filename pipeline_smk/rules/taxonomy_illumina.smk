import os
rule illumina_kraken2:
    input:
        dbjson=DB_JSON,
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R2.fastq.gz")
    output:
        kraken_out=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "kraken2_output.txt"),
        kraken_report=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "kraken2_report.txt")
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_kraken2.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_kraken2.txt")
    conda:
        os.path.join(ENVS_DIR, "kraken.yaml")
    threads: THREADS
    resources:
        mem_mb=16000
    params:
        confidence=lambda wc: config["taxonomy"]["kraken2_confidence"],
        mem_map=lambda wc: "--memory_mapping" if config["taxonomy"].get("kraken2_memory_mapping", False) else ""
    shell:
        r"""
        exec > {log} 2>&1
        python3 scripts/run_kraken2.py \
          --dbjson {input.dbjson} \
          --out {output.kraken_out} \
          --report {output.kraken_report} \
          --threads {threads} \
          --confidence {params.confidence} \
          {params.mem_map} \
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
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_tb_flag.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_tb_flag.txt")
    shell:
        r"""
        exec > {log} 2>&1
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
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_tbprofiler.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_tbprofiler.txt")
    conda:
        os.path.join(ENVS_DIR, "tbprofiler.yaml")
    threads: THREADS
    shell:
        r"""
        exec > {log} 2>&1
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
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_centrifuge.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_centrifuge.txt")
    conda:
        os.path.join(ENVS_DIR, "centrifuge.yaml")
    threads: THREADS
    params:
        min_hitlen=lambda wc: config["taxonomy"]["centrifuge_min_hitlen_illumina"],
        min_score=lambda wc: config["taxonomy"].get("centrifuge_min_score", 0)
    shell:
        r"""
        exec > {log} 2>&1
        python3 scripts/run_centrifuge.py \
          --dbjson {input.dbjson} \
          --out {output.centrifuge_out} \
          --report {output.centrifuge_report} \
          --kreport {output.centrifuge_kreport} \
          --threads {threads} \
          --min_hitlen {params.min_hitlen} \
          --min_score {params.min_score} \
          --r1 {input.r1} --r2 {input.r2}
        """
rule illumina_krona_kraken2:
    input:
        kraken_report=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "kraken2_report.txt"),
        dbjson=DB_JSON
    output:
        html=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "krona_kraken2.html")
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_krona_kraken2.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_krona_kraken2.txt")
    conda:
        os.path.join(ENVS_DIR, "krona.yaml")
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p $(dirname {output.html})
        python3 scripts/kreport2text.py {input.kraken_report} {output.html}.txt
        ktImportText -o {output.html} {output.html}.txt
        rm -f {output.html}.txt
        """
rule illumina_krona_centrifuge:
    input:
        centrifuge_kreport=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_kreport.txt"),
        dbjson=DB_JSON
    output:
        html=os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "krona_centrifuge.html")
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_krona_centrifuge.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_krona_centrifuge.txt")
    conda:
        os.path.join(ENVS_DIR, "krona.yaml")
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p $(dirname {output.html})
        python3 scripts/kreport2text.py {input.centrifuge_kreport} {output.html}.txt
        ktImportText -o {output.html} {output.html}.txt
        rm -f {output.html}.txt
        """
