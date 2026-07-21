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

rule nanopore_kma_resfinder:
    input:
        r1=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "nonhost_nanopore.fastq.gz"),
        amr_mode=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "amr_mode.txt"),
        dbjson=DB_JSON
    output:
        kma_res=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "kma_resfinder", "kma.res"),
        kma_tsv=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "kma_resfinder", "kma.tsv")
    conda:
        os.path.join(ENVS_DIR, "kma.yaml")
    threads: THREADS
    shell:
        r"""
        mkdir -p $(dirname {output.kma_res})
        python3 scripts/run_kma_resfinder.py \
          --dbjson {input.dbjson} \
          --amr_mode_file {input.amr_mode} \
          --r1 {input.r1} \
          --platform nanopore \
          --outdir $(dirname {output.kma_res}) --threads {threads}
        """
