import os

rule nanopore_contig_amr:
    input:
        contigs = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contigs.fasta")
    output:
        abricate_card = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "abricate_card.tsv"),
        amrfinder = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "amrfinder_results.tsv")
    conda:
        os.path.join(ENVS_DIR, "amr_contigs.yaml")
    threads: THREADS
    shell:
        r"""
        python3 scripts/run_contig_amr.py \
          --contigs {input.contigs} \
          --threads {threads} \
          --minid {config[amr][abricate_min_identity]} \
          --mincov {config[amr][abricate_min_coverage]} \
          --outdir $(dirname {output.abricate_card})
        """
