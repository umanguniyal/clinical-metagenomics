import os
rule run_prokka_bins:
    input:
        bin_dir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning"),
        checkm2_report = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_checkm2", "quality_report.tsv")
    output:
        out_dir = directory(os.path.join(RESULTS_DIR, "{platform}", "{sample}", "bin_annotations")),
        stats = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "bin_annotations", "prokka_stats.json")
    log:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "logs", "run_prokka_bins.log")
    benchmark:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "benchmarks", "run_prokka_bins.txt")
    conda:
        os.path.join(ENVS_DIR, "prokka.yaml")
    threads: THREADS
    shell:
        r"""
        exec > {log} 2>&1
        python3 scripts/run_prokka_bins.py \
          --bin_dir {input.bin_dir} \
          --checkm2_report {input.checkm2_report} \
          --out_dir {output.out_dir} \
          --threads {threads}
        """
