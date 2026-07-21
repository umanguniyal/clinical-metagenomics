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
        heatmap     = os.path.join(RESULTS_DIR, "illumina", "{sample}", "ml_results", "amr_heatmap.png"),
        kma_scatter = os.path.join(RESULTS_DIR, "illumina", "{sample}", "ml_results", "kma_precision_scatter.png"),
        score_bar   = os.path.join(RESULTS_DIR, "illumina", "{sample}", "ml_results", "pathogen_scores_bar.png"),
        pdf         = os.path.join(RESULTS_DIR, "illumina", "{sample}", "report", "report_bioinformatician.pdf")
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
        heatmap     = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "ml_results", "amr_heatmap.png"),
        kma_scatter = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "ml_results", "kma_precision_scatter.png"),
        score_bar   = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "ml_results", "pathogen_scores_bar.png"),
        pdf         = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "report", "report_bioinformatician.pdf")
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
        python3 scripts/generate_pdf_report.py \
          --report_json {input.report_json} \
          --out         {output.pdf}
        """

rule cohort_ml_analysis:
    input:
        reports = illumina_reports() + nanopore_reports()
    output:
        pca = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_pca_landscape.png"),
        burden = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_amr_burden.png"),
        matrix = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_species_amr_matrix.png"),
        dist = os.path.join(RESULTS_DIR, "cohort_ml_results", "cohort_score_distribution.png")
    conda:
        os.path.join(ENVS_DIR, "ml_metagenomics.yaml")
    threads: 2
    shell:
        r"""
        mkdir -p $(dirname {output.pca})
        if [ "{input.reports}" != "" ]; then
            python3 scripts/aggregate_cohort_ml.py \
                --reports {input.reports} \
                --outdir $(dirname {output.pca})
        else
            touch {output.pca} {output.burden} {output.matrix} {output.dist}
        fi
        """
