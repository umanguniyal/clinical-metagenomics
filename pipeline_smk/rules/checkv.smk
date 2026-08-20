import os
rule run_checkv:
    input:
        dbjson = DB_JSON,
        bins_dir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_viral")
    output:
        quality = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_checkv", "quality_summary.tsv")
    log:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "logs", "run_checkv.log")
    benchmark:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "benchmarks", "run_checkv.txt")
    conda:
        os.path.join(ENVS_DIR, "checkv.yaml")
    threads: THREADS
    resources:
        mem_mb=8000
    params:
        checkv_outdir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_checkv")
    shell:
        r"""
        exec > {log} 2>&1
        python3 scripts/run_checkv.py --dbjson {input.dbjson} --bin_dir {input.bins_dir} --threads {threads} --outdir {params.checkv_outdir}
        """
