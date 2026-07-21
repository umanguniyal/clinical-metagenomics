import os

# ==== REPORT GENERATION RULES ====

rule nanopore_clinician_report:
    input:
        mean_cov = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "mean_coverage.txt"),
        amr_mode = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "amr_mode.txt"),
        kraken_report = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "kraken2_report.txt"),
        centrifuge_out = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "centrifuge_output.txt"),
        centrifuge_report = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "centrifuge_report.tsv"),
        centrifuge_kreport = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "centrifuge_kreport.txt"),
        quast_tsv = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "quast", "report.tsv"),
        abricate_card = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "abricate_card.tsv"),
        amrfinder = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "amrfinder_results.tsv"),
        plasmid_fasta = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "plasmids", "plasmid_unbinned.fasta"),
        mobtyper = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "plasmids", "mobtyper_results.txt"),
        contig_report = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "plasmids", "contig_report.txt"),
        bins_dir = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "binning"),
        kma_res = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "kma_resfinder", "kma.res"),
        kma_tsv = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "kma_resfinder", "kma.tsv"),
        stats = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "host_stats.json")
    output:
        json = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "report", "report_clinician.json"),
        txt = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "report", "report_clinician.txt"),
        pdf = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "report", "report_clinician.pdf")
    conda:
        os.path.join(ENVS_DIR, "base_utils.yaml")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.json})
        python3 scripts/make_report_clinician.py \
          --platform nanopore \
          --profile clinician \
          --sample {wildcards.sample} \
          --mean_cov {input.mean_cov} \
          --amr_mode {input.amr_mode} \
          --kraken_report {input.kraken_report} \
          --bracken_report "" \
          --centrifuge_out {input.centrifuge_out} \
          --centrifuge_report {input.centrifuge_report} \
          --centrifuge_kreport {input.centrifuge_kreport} \
          --quast_tsv {input.quast_tsv} \
          --abricate_card {input.abricate_card} \
          --amrfinder {input.amrfinder} \
          --plasmid_fasta {input.plasmid_fasta} \
          --mobtyper_results {input.mobtyper} \
          --contig_report {input.contig_report} \
          --bins_dir {input.bins_dir} \
          --kma_res {input.kma_res} \
          --kma_tsv {input.kma_tsv} \
          --host_stats {input.stats} \
          --out {output.json} \
          --out_txt {output.txt}
        python3 scripts/generate_pdf_report.py \
          --report_json {output.json} \
          --out {output.pdf}
        """

rule nanopore_bioinformatician_report:
    input:
        mean_cov = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "mean_coverage.txt"),
        amr_mode = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "amr_mode.txt"),
        kraken_report = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "kraken2_report.txt"),
        centrifuge_out = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "centrifuge_output.txt"),
        centrifuge_report = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "centrifuge_report.tsv"),
        centrifuge_kreport = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "taxonomy", "centrifuge_kreport.txt"),
        quast_tsv = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "quast", "report.tsv"),
        abricate_card = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "abricate_card.tsv"),
        amrfinder = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "amrfinder_results.tsv"),
        bin_amr_json = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "bin_amr_map.json"),
        plasmid_fasta = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "plasmids", "plasmid_unbinned.fasta"),
        mobtyper = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "plasmids", "mobtyper_results.txt"),
        contig_report = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "plasmids", "contig_report.txt"),
        bins_dir = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "binning"),
        kma_res = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "kma_resfinder", "kma.res"),
        kma_tsv = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "amr", "kma_resfinder", "kma.tsv"),
        stats = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "host_stats.json")
    output:
        json = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "report", "report_bioinformatician.json"),
        txt = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "report", "report_bioinformatician.txt")
    conda:
        os.path.join(ENVS_DIR, "base_utils.yaml")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.json})
        python3 scripts/make_report.py \
          --platform nanopore \
          --profile bioinformatician \
          --sample {wildcards.sample} \
          --mean_cov {input.mean_cov} \
          --amr_mode {input.amr_mode} \
          --kraken_report {input.kraken_report} \
          --bracken_report "" \
          --centrifuge_out {input.centrifuge_out} \
          --centrifuge_report {input.centrifuge_report} \
          --centrifuge_kreport {input.centrifuge_kreport} \
          --quast_tsv {input.quast_tsv} \
          --abricate_card {input.abricate_card} \
          --amrfinder {input.amrfinder} \
          --bin_amr_json {input.bin_amr_json} \
          --plasmid_fasta {input.plasmid_fasta} \
          --mobtyper_results {input.mobtyper} \
          --contig_report {input.contig_report} \
          --bins_dir {input.bins_dir} \
          --kma_res {input.kma_res} \
          --kma_tsv {input.kma_tsv} \
          --host_stats {input.stats} \
          --out {output.json} \
          --out_txt {output.txt}
        """
