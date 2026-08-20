import os
rule nanopore_assemble_flye:
    input:
        fq = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "host_removed", "nonhost_nanopore.fastq.gz"),
        status=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "sample_status.json")
    output:
        contigs = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contigs.fasta")
    log:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "logs", "nanopore_assemble_flye.log")
    benchmark:
        os.path.join(RESULTS_DIR, "nanopore", "{sample}", "benchmarks", "nanopore_assemble_flye.txt")
    conda:
        os.path.join(ENVS_DIR, "flye.yaml")
    threads: THREADS
    resources:
        mem_mb=16000
    params:
        assembly_dir=os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly")
    shell:
        r"""
        exec > {log} 2>&1
        mkdir -p {params.assembly_dir}
        
        sufficient=$(python3 -c "import json; print(json.load(open('{input.status}'))['sufficient'])" 2>/dev/null || echo "False")
        if [ "$sufficient" == "False" ]; then
            echo "Skipping assembly due to insufficient reads/bases"
            touch {output.contigs}
            exit 0
        fi
        
        FLYE_TMP={params.assembly_dir}/.flye_tmp_$$
        rm -rf $FLYE_TMP
        ( flye --nano-raw {input.fq} --out-dir $FLYE_TMP \
             --meta \
             --min-overlap {config[assembly][flye_min_overlap]} \
             --threads {threads} ) || true
             
        if [ ! -f "$FLYE_TMP/assembly.fasta" ] || [ ! -s "$FLYE_TMP/assembly.fasta" ]; then
            echo "WARNING: Flye produced no contigs for {wildcards.sample}" >&2
            echo -e "{wildcards.sample}\tassembly\tFlye produced no contigs or failed" >> {params.assembly_dir}/../../../pipeline_errors.log
            touch {output.contigs}
        else
            cp $FLYE_TMP/assembly.fasta {output.contigs}
        fi
        rm -rf $FLYE_TMP
        """
