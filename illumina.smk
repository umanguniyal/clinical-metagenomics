import os

# ==== REPORT GENERATION RULES ====

rule illumina_clinician_report:
    input:
        mean_cov = os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "mean_coverage.txt"),
        amr_mode = os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "amr_mode.txt"),
        kraken_report = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "kraken2_report.txt"),
        bracken_report = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "bracken_report.tsv"),
        centrifuge_out = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_output.txt"),
        centrifuge_report = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_report.tsv"),
        centrifuge_kreport = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_kreport.txt"),
        quast_tsv = os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "quast", "report.tsv"),
        abricate_card = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "abricate_card.tsv"),
        amrfinder = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "amrfinder_results.tsv"),
        kma_res = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "kma_resfinder", "kma.res"),
        kma_tsv = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "kma_resfinder", "kma.tsv")
    output:
        json = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_clinician.json"),
        txt = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_clinician.txt")
    conda:
        os.path.join(ENVS_DIR, "base_utils.yaml")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.json})
        python3 scripts/make_report_clinician.py \
          --platform illumina \
          --profile clinician \
          --sample {wildcards.sample} \
          --mean_cov {input.mean_cov} \
          --amr_mode {input.amr_mode} \
          --kraken_report {input.kraken_report} \
          --bracken_report {input.bracken_report} \
          --centrifuge_out {input.centrifuge_out} \
          --centrifuge_report {input.centrifuge_report} \
          --centrifuge_kreport {input.centrifuge_kreport} \
          --quast_tsv {input.quast_tsv} \
          --abricate_card {input.abricate_card} \
          --amrfinder {input.amrfinder} \
          --kma_res {input.kma_res} \
          --kma_tsv {input.kma_tsv} \
          --out {output.json} \
          --out_txt {output.txt}
        """

rule illumina_bioinformatician_report:
    input:
        mean_cov = os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "mean_coverage.txt"),
        amr_mode = os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "amr_mode.txt"),
        kraken_report = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "kraken2_report.txt"),
        bracken_report = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "bracken_report.tsv"),
        centrifuge_out = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_output.txt"),
        centrifuge_report = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_report.tsv"),
        centrifuge_kreport = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "centrifuge_kreport.txt"),
        quast_tsv = os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "quast", "report.tsv"),
        abricate_card = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "abricate_card.tsv"),
        amrfinder = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "amrfinder_results.tsv"),
        bin_amr_json = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "bin_amr_map.json")
    output:
        json = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_bioinformatician.json"),
        txt = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_bioinformatician.txt")
    conda:
        os.path.join(ENVS_DIR, "base_utils.yaml")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.json})
        python3 scripts/make_report.py \
          --platform illumina \
          --profile bioinformatician \
          --sample {wildcards.sample} \
          --mean_cov {input.mean_cov} \
          --amr_mode {input.amr_mode} \
          --kraken_report {input.kraken_report} \
          --bracken_report {input.bracken_report} \
          --centrifuge_out {input.centrifuge_out} \
          --centrifuge_report {input.centrifuge_report} \
          --centrifuge_kreport {input.centrifuge_kreport} \
          --quast_tsv {input.quast_tsv} \
          --abricate_card {input.abricate_card} \
          --amrfinder {input.amrfinder} \
          --bin_amr_json {input.bin_amr_json} \
          --out {output.json} \
          --out_txt {output.txt}
        """
