import os

def get_control_json_illumina(wildcards):
    if wildcards.sample + "_control" in SAMPLES_ILLUMINA:
        return os.path.join(RESULTS_DIR, "illumina", wildcards.sample + "_control", "report", f"report_{PROFILE}.json")
    return []

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
        kma_tsv = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "kma_resfinder", "kma.tsv"),
        plasmid_fasta = os.path.join(RESULTS_DIR, "illumina", "{sample}", "plasmids", "plasmid_unbinned.fasta"),
        mobtyper = os.path.join(RESULTS_DIR, "illumina", "{sample}", "plasmids", "mobtyper_results.txt"),
        contig_report = os.path.join(RESULTS_DIR, "illumina", "{sample}", "plasmids", "contig_report.txt"),
        bins_dir = os.path.join(RESULTS_DIR, "illumina", "{sample}", "binning"),
        stats = os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "host_stats.json"),
        tbprofiler_json = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "tbprofiler", "results", "{sample}.results.json"),
        control_json = get_control_json_illumina
    output:
        json = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_clinician.json"),
        txt = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_clinician.txt"),
        pdf = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_clinician.pdf")
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_clinician_report.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_clinician_report.txt")
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
          --plasmid_fasta {input.plasmid_fasta} \
          --mobtyper_results {input.mobtyper} \
          --contig_report {input.contig_report} \
          --bins_dir {input.bins_dir} \
          --host_stats {input.stats} \
          --tbprofiler_json {input.tbprofiler_json} \
          --checkv_tsv {input.checkv_tsv} \
          --busco_txt {input.busco_txt} \
          --out {output.json} \
          --out_txt {output.txt}
        
        if [ -n "{input.control_json}" ] && [ -s "{input.control_json}" ]; then
            python3 scripts/qc_crosscheck.py --patient_json {output.json} --control_json {input.control_json} --out_json {output.json}.tmp
            mv {output.json}.tmp {output.json}
        fi
        
        python3 scripts/generate_pdf_report.py \
          --report_json {output.json} \
          --out {output.pdf}
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
        bin_amr_json = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "bin_amr_map.json"),
        plasmid_fasta = os.path.join(RESULTS_DIR, "illumina", "{sample}", "plasmids", "plasmid_unbinned.fasta"),
        mobtyper = os.path.join(RESULTS_DIR, "illumina", "{sample}", "plasmids", "mobtyper_results.txt"),
        contig_report = os.path.join(RESULTS_DIR, "illumina", "{sample}", "plasmids", "contig_report.txt"),
        bins_dir = os.path.join(RESULTS_DIR, "illumina", "{sample}", "binning"),
        kma_res = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "kma_resfinder", "kma.res"),
        kma_tsv = os.path.join(RESULTS_DIR, "illumina", "{sample}", "amr", "kma_resfinder", "kma.tsv"),
        prokka_json = os.path.join(RESULTS_DIR, "illumina", "{sample}", "bin_annotations", "prokka_stats.json"),
        stats = os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "host_stats.json"),
        tbprofiler_json = os.path.join(RESULTS_DIR, "illumina", "{sample}", "taxonomy", "tbprofiler", "results", "{sample}.results.json"),
        checkv_tsv = os.path.join(RESULTS_DIR, "illumina", "{sample}", "binning_checkv", "quality_summary.tsv"),
        busco_txt = os.path.join(RESULTS_DIR, "illumina", "{sample}", "binning_busco", "short_summary.txt"),
        control_json = get_control_json_illumina
    output:
        json = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_bioinformatician.json"),
        txt = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_bioinformatician.txt")
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_bioinformatician_report.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_bioinformatician_report.txt")
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
          --plasmid_fasta {input.plasmid_fasta} \
          --mobtyper_results {input.mobtyper} \
          --contig_report {input.contig_report} \
          --bins_dir {input.bins_dir} \
          --kma_res {input.kma_res} \
          --kma_tsv {input.kma_tsv} \
          --prokka_json {input.prokka_json} \
          --host_stats {input.stats} \
          --tbprofiler_json {input.tbprofiler_json} \
          --checkv_tsv {input.checkv_tsv} \
          --busco_txt {input.busco_txt} \
          --out {output.json} \
          --out_txt {output.txt}
        
        if [ -n "{input.control_json}" ] && [ -s "{input.control_json}" ]; then
            python3 scripts/qc_crosscheck.py --patient_json {output.json} --control_json {input.control_json} --out_json {output.json}.tmp
            mv {output.json}.tmp {output.json}
        fi
        """
