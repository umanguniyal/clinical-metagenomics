import os

rule illumina_assemble_megahit:
    input:
        r1=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R1.fastq.gz"),
        r2=os.path.join(RESULTS_DIR, "illumina", "{sample}", "host_removed", "nonhost_R2.fastq.gz"),
        status=os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "sample_status.json")
    output:
        contigs=os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "contigs.fasta")
    log:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "logs", "illumina_assemble_megahit.log")
    benchmark:
        os.path.join(RESULTS_DIR, "illumina", "{sample}", "benchmarks", "illumina_assemble_megahit.txt")
    conda:
        os.path.join(ENVS_DIR, "megahit.yaml")
    threads: THREADS
    resources:
        mem_mb=16000
    params:
        assembly_dir=os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly")
    shell:
        r"""
        mkdir -p {params.assembly_dir}
        
        sufficient=$(python3 -c "import json; print(json.load(open('{input.status}'))['sufficient'])" 2>/dev/null || echo "False")
        if [ "$sufficient" == "False" ]; then
            echo "Skipping assembly due to insufficient reads/bases"
            touch {output.contigs}
            exit 0
        fi
        
        MEGAHIT_TMP={params.assembly_dir}/.megahit_tmp_$$
        rm -rf $MEGAHIT_TMP
        ( megahit -1 {input.r1} -2 {input.r2} \
            -o $MEGAHIT_TMP \
            --min-contig-len {config[assembly][megahit_min_contig_len]} \
            --k-min {config[assembly][megahit_k_min]} \
            --k-max {config[assembly][megahit_k_max]} \
            --k-step {config[assembly][megahit_k_step]} \
            --num-cpu-threads {threads} ) || true
            
        if [ ! -f "$MEGAHIT_TMP/final.contigs.fa" ] || [ ! -s "$MEGAHIT_TMP/final.contigs.fa" ]; then
            echo "WARNING: MEGAHIT produced no contigs for {wildcards.sample}" >&2
            echo -e "{wildcards.sample}\tassembly\tMEGAHIT produced no contigs or failed" >> {params.assembly_dir}/../../../pipeline_errors.log
            touch {output.contigs}
        else
            cp $MEGAHIT_TMP/final.contigs.fa {output.contigs}
        fi
        rm -rf $MEGAHIT_TMP
        """
