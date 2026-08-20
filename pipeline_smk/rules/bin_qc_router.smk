import os
rule bin_qc_router:
    input:
        bins_dir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning"),
        tax_dir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "bin_taxonomy")
    output:
        bact_dir = directory(os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_bacteria")),
        viral_dir = directory(os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_viral")),
        euk_dir = directory(os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning_eukarya"))
    log:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "logs", "bin_qc_router.log")
    benchmark:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "benchmarks", "bin_qc_router.txt")
    threads: 1
    shell:
        r"""
        exec > {log} 2>&1
        python3 scripts/bin_qc_router.py \
            --bin_dir {input.bins_dir} \
            --bin_tax_dir {input.tax_dir} \
            --out_checkm2 {output.bact_dir} \
            --out_checkv {output.viral_dir} \
            --out_busco {output.euk_dir}
        """
