import os

rule illumina_assemble_megahit:
    input:
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R2.fastq.gz")
    output:
        contigs=os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "contigs.fasta")
    conda:
        os.path.join(ENVS_DIR, "megahit.yaml")
    threads: THREADS
    params:
        assembly_dir=os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly")
    shell:
        r"""
        mkdir -p {params.assembly_dir}
        MEGAHIT_TMP={params.assembly_dir}/.megahit_tmp_$$
        rm -rf $MEGAHIT_TMP
        megahit -1 {input.r1} -2 {input.r2} \
            -o $MEGAHIT_TMP \
            --min-contig-len {config[assembly][megahit_min_contig_len]} \
            --k-min {config[assembly][megahit_k_min]} \
            --k-max {config[assembly][megahit_k_max]} \
            --k-step {config[assembly][megahit_k_step]} \
            --num-cpu-threads {threads}
        cp $MEGAHIT_TMP/final.contigs.fa {output.contigs}
        rm -rf $MEGAHIT_TMP
        """
