import os

rule run_checkm2:
    input:
        dbjson = DB_JSON,
        bins_dir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning")
    output:
        quality = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_checkm2", "quality_report.tsv")
    conda:
        os.path.join(ENVS_DIR, "checkm2.yaml")
    threads: THREADS
    params:
        checkm2_outdir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_checkm2")
    shell:
        r"""
        python3 scripts/run_checkm2.py --dbjson {input.dbjson} --bin_dir {input.bins_dir} --threads {threads} --outdir {params.checkm2_outdir}
        """
