import os
rule run_checkm2:
    input:
        dbjson = DB_JSON,
        bins_dir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_bacteria")
    output:
        quality = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_checkm2", "quality_report.tsv")
    log:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "logs", "run_checkm2.log")
    benchmark:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "benchmarks", "run_checkm2.txt")
    conda:
        os.path.join(ENVS_DIR, "checkm2.yaml")
    threads: THREADS
    resources:
        mem_mb=16000
    params:
        checkm2_outdir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_checkm2")
    shell:
        r"""
        exec > {log} 2>&1
        python3 scripts/run_checkm2.py --dbjson {input.dbjson} --bin_dir {input.bins_dir} --threads {threads} --outdir {params.checkm2_outdir}
        """
