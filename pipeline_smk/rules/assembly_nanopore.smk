import os

rule nanopore_assemble_flye:
    input:
        fq = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "nonhost_nanopore.fastq.gz")
    output:
        contigs = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contigs.fasta")
    conda:
        os.path.join(ENVS_DIR, "flye.yaml")
    threads: THREADS
    params:
        assembly_dir=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly")
    shell:
        r"""
        mkdir -p {params.assembly_dir}
        FLYE_TMP={params.assembly_dir}/.flye_tmp_$$
        rm -rf $FLYE_TMP
        flye --nano-raw {input.fq} --out-dir $FLYE_TMP \
             --min-overlap {config[assembly][flye_min_overlap]} \
             --threads {threads} \
             --genome-size 4m
        cp $FLYE_TMP/assembly.fasta {output.contigs}
        rm -rf $FLYE_TMP
        """
