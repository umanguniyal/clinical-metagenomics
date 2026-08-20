import os

rule run_busco:
    input:
        dbjson = DB_JSON,
        bins_dir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_eukarya")
    output:
        quality = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_busco", "short_summary.txt")
    log:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "logs", "run_busco.log")
    benchmark:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "benchmarks", "run_busco.txt")
    conda:
        os.path.join(ENVS_DIR, "busco.yaml")
    threads: THREADS
    resources:
        mem_mb=8000
    params:
        busco_outdir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_busco")
    shell:
        r"""
        python3 scripts/run_busco.py --dbjson {input.dbjson} --bin_dir {input.bins_dir} --threads {threads} --outdir {params.busco_outdir}
        """
