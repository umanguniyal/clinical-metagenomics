import os

# ── ML Analysis rules — bioinformatician profile ONLY ──────────────────────────
# These rules are only included when PROFILE == "bioinformatician" (see Snakefile)

rule illumina_ml_analysis:
    input:
        report_json = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_bioinformatician.json")
    output:
        predictions = os.path.join(RESULTS_DIR, "illumina", "{sample}", "ml_results", "predictions.csv"),
        feat_imp    = os.path.join(RESULTS_DIR, "illumina", "{sample}", "ml_results", "feature_importance.png"),
        pca         = os.path.join(RESULTS_DIR, "illumina", "{sample}", "ml_results", "pca_clusters.png"),
        heatmap     = os.path.join(RESULTS_DIR, "illumina", "{sample}", "ml_results", "organism_quality_heatmap.png"),
        kma_scatter = os.path.join(RESULTS_DIR, "illumina", "{sample}", "ml_results", "kma_precision_scatter.png"),
        score_bar   = os.path.join(RESULTS_DIR, "illumina", "{sample}", "ml_results", "pathogen_scores_bar.png"),
        checkm2     = os.path.join(RESULTS_DIR, "illumina", "{sample}", "ml_results", "checkm2_vs_prokka.png"),
        pdf         = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_bioinformatician.pdf")
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_ml_analysis.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_ml_analysis.txt")
    conda:
        os.path.join(ENVS_DIR, "ml_metagenomics.yaml")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.predictions})
        python3 scripts/ml_pathogen_amr.py \
          --report_json {input.report_json} \
          --outdir      $(dirname {output.predictions}) \
          --model       rf
        python3 scripts/bioinfo_triage_lists.py \
          --report_json {input.report_json} \
          --bartlett_xlsx bacteria_human_pathogens.xlsx \
          --viral_xlsx viral_pathogens.xlsx \
          --fungal_xlsx fungal_pathogens.xlsx \
          --out_txt $(dirname {input.report_json})/report_bioinformatician.txt
        python3 scripts/generate_pdf_report.py \
          --report_json {input.report_json} \
          --out         {output.pdf}
        """

rule nanopore_ml_analysis:
    input:
        report_json = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "report", "report_bioinformatician.json")
    output:
        predictions = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "ml_results", "predictions.csv"),
        feat_imp    = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "ml_results", "feature_importance.png"),
        pca         = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "ml_results", "pca_clusters.png"),
        heatmap     = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "ml_results", "organism_quality_heatmap.png"),
        kma_scatter = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "ml_results", "kma_precision_scatter.png"),
        score_bar   = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "ml_results", "pathogen_scores_bar.png"),
        checkm2     = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "ml_results", "checkm2_vs_prokka.png"),
        pdf         = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "report", "report_bioinformatician.pdf")
    log:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "logs", "nanopore_ml_analysis.log")
    benchmark:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "benchmarks", "nanopore_ml_analysis.txt")
    conda:
        os.path.join(ENVS_DIR, "ml_metagenomics.yaml")
    threads: 4
    shell:
        r"""
        mkdir -p $(dirname {output.predictions})
        python3 scripts/ml_pathogen_amr.py \
          --report_json {input.report_json} \
          --outdir      $(dirname {output.predictions}) \
          --model       rf
        python3 scripts/bioinfo_triage_lists.py \
          --report_json {input.report_json} \
          --bartlett_xlsx bacteria_human_pathogens.xlsx \
          --viral_xlsx viral_pathogens.xlsx \
          --fungal_xlsx fungal_pathogens.xlsx \
          --out_txt $(dirname {input.report_json})/report_bioinformatician.txt
        python3 scripts/generate_pdf_report.py \
          --report_json {input.report_json} \
          --out         {output.pdf}
        """

rule cohort_ml_analysis:
    input:
        reports = glob.glob(os.path.join(RESULTS_DIR, "*", "*", "report", "report_bioinformatician.json"))
    output:
        pca = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_pca_landscape.png"),
        burden = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_amr_burden.png"),
        matrix = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_species_amr_matrix.png"),
        dist = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_score_distribution.png"),
        c_checkm2 = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_checkm2_vs_prokka.png"),
        c_stab = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_feature_stability.png"),
        c_cal = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_calibration_curve.png")
    log:
        os.path.join(RESULTS_DIR, "logs", "cohort_ml_analysis.log")
    benchmark:
        os.path.join(RESULTS_DIR, "benchmarks", "cohort_ml_analysis.txt")
    conda:
        os.path.join(ENVS_DIR, "ml_metagenomics.yaml")
    threads: 2
    shell:
        r"""
        mkdir -p $(dirname {output.pca})
        count=$(echo "{input.reports}" | wc -w)
        if [ "$count" -ge 1 ]; then
            python3 scripts/aggregate_cohort_ml.py \
                --reports {input.reports} \
                --outdir $(dirname {output.pca})
        else
            echo "No completed reports found. Touching dummy outputs."
            touch {output.pca} {output.burden} {output.matrix} {output.dist} {output.c_checkm2} {output.c_stab} {output.c_cal}
        fi
        """
