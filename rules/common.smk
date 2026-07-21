import os

def _sample_fastq_for_readlen(wc):
    if wc.platform == "illumina":
        return ILLUMINA_SAMPLES[wc.sample]["r1"]
    return NANOPORE_SAMPLES[wc.sample]["fq"]

rule estimate_read_length:
    input:
        fastq=_sample_fastq_for_readlen,
        dbjson=DB_JSON
    output:
        os.path.join(RESULTS_DIR, "{platform}", "{sample}", "qc", "estimated_read_length.txt")
    conda:
        os.path.join(ENVS_DIR, "base_utils.yaml")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output})
        python3 scripts/estimate_read_length.py --fastq {input.fastq} --out {output}
        """

# ONLY RUN BRACKEN FOR ILLUMINA
if PLATFORM == "illumina":
    rule bracken:
        input:
            kraken_report=os.path.join(RESULTS_DIR, "{platform}", "{sample}", "taxonomy", "kraken2_report.txt"),
            readlen=os.path.join(RESULTS_DIR, "{platform}", "{sample}", "qc", "estimated_read_length.txt"),
            dbjson=DB_JSON
        output:
            bracken_report=os.path.join(RESULTS_DIR, "{platform}", "{sample}", "taxonomy", "bracken_report.tsv"),
            bracken_species=os.path.join(RESULTS_DIR, "{platform}", "{sample}", "taxonomy", "bracken_species.tsv")
        conda:
            os.path.join(ENVS_DIR, "kraken.yaml")
        threads: THREADS
        shell:
            r"""
            mkdir -p $(dirname {output.bracken_report})
            python3 scripts/run_bracken.py \
              --dbjson {input.dbjson} \
              --platform {wildcards.platform} \
              --readlen_file {input.readlen} \
              --kraken_report {input.kraken_report} \
              --out_report {output.bracken_report} \
              --out_species {output.bracken_species} \
              --level "{config[taxonomy][bracken_level]}"
            """
