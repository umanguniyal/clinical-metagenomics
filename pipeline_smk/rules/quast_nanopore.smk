import os

rule nanopore_quast:
    input:
        contigs = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contigs.fasta")
    output:
        report = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "quast", "report.tsv")
    conda:
        os.path.join(ENVS_DIR, "quast.yaml")
    threads: THREADS
    params:
        quast_dir=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "quast")
    shell:
        r"""
        mkdir -p {params.quast_dir}
        rm -rf {params.quast_dir}/.quast_tmp
        quast.py {input.contigs} -o {params.quast_dir}/.quast_tmp --threads {threads}
        cp {params.quast_dir}/.quast_tmp/report.tsv {output.report}
        rm -rf {params.quast_dir}/.quast_tmp
        """
