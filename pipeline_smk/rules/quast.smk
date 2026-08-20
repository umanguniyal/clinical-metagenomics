import os
rule illumina_quast:
    input:
        contigs=os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "contigs.fasta")
    output:
        report=os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "quast", "report.tsv")
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_quast.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_quast.txt")
    conda:
        os.path.join(ENVS_DIR, "quast.yaml")
    threads: THREADS
    params:
        quast_dir=os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "quast")
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.quast_dir}
        rm -rf {params.quast_dir}/.quast_tmp
        
        if [ ! -s {input.contigs} ] || ! grep -q "^>" {input.contigs}; then
            echo -e "Assembly\tcontigs\n# contigs\t0\nTotal length\t0\nLargest contig\t0\nN50\t0\nGC (%)\t0" > {output.report}
        else
            quast.py {input.contigs} -o {params.quast_dir}/.quast_tmp --threads {threads} || true
            if [ -f {params.quast_dir}/.quast_tmp/report.tsv ]; then
                cp {params.quast_dir}/.quast_tmp/report.tsv {output.report}
            else
                echo -e "Assembly\tcontigs\n# contigs\t0\nTotal length\t0\nLargest contig\t0\nN50\t0\nGC (%)\t0" > {output.report}
            fi
            rm -rf {params.quast_dir}/.quast_tmp
        fi
        """
