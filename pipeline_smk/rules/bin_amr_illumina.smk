import os

rule illumina_bin_amr_map:
    input:
        bin_dir=os.path.join(RESULTS_DIR, "illumina", "{sample}", "binning"),
        bin_tax_dir=os.path.join(RESULTS_DIR, "illumina", "{sample}", "bin_taxonomy"),
        contig_abricate_card=os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "abricate_card.tsv"),
        contig_amrfinder=os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "amrfinder_results.tsv"),
        checkm2_quality=os.path.join(RESULTS_DIR, "illumina", "{sample}", "binning_checkm2", "quality_report.tsv")
    output:
        out_json=os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "bin_amr_map.json"),
        out_tsv=os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "bin_amr_map.tsv"),
        out_summary=os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "bin_amr_map_summary.tsv")
    conda:
        os.path.join(ENVS_DIR, "bin_amr.yaml")
    threads: THREADS
    shell:
        r"""
        python3 scripts/bin_amr_map.py \
          --bin_dir {input.bin_dir} \
          --bin_taxonomy_dir {input.bin_tax_dir} \
          --contig_abricate_card {input.contig_abricate_card} \
          --contig_amrfinder {input.contig_amrfinder} \
          --checkm2_quality {input.checkm2_quality} \
          --out_json {output.out_json} --out_tsv {output.out_tsv} --out_summary {output.out_summary}
        """
