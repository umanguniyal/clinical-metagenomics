import os

rule run_prokka_bins:
    input:
        bin_dir = os.path.join(RESULTS_DIR, "{platform}", "{sample}", "binning")
    output:
        out_dir = directory(os.path.join(RESULTS_DIR, "{platform}", "{sample}", "bin_annotations"))
    conda:
        os.path.join(ENVS_DIR, "prokka.yaml")
    threads: THREADS
    shell:
        r"""
        python3 scripts/run_prokka_bins.py \
          --bin_dir {input.bin_dir} --out_dir {output.out_dir} --threads {threads}
        """ 
