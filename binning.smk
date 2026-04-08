import os

rule illumina_binning_bam:
    input:
        contigs = os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "contigs.fasta"),
        reads1 = lambda wc: ILLUMINA_SAMPLES[wc.sample]["r1"],
        reads2 = lambda wc: ILLUMINA_SAMPLES[wc.sample]["r2"]
    output:
        bam = os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "contigs_sorted.bam"),
        bai = os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "contigs_sorted.bam.bai")
    conda:
        os.path.join(ENVS_DIR, "coverage_illumina.yaml")
    threads: THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        
        # Build bowtie2 index
        TMPIDX=$(mktemp -d)
        bowtie2-build {input.contigs} $TMPIDX/idx --threads {threads} 2>/dev/null
        
        # Map reads and create sorted BAM
        bowtie2 -x $TMPIDX/idx -1 {input.reads1} -2 {input.reads2} --threads {threads} -q 2>/dev/null | \
          samtools view -bS -q 30 - | \
          samtools sort -@ {threads} -o {output.bam}
        
        # Index BAM file
        samtools index {output.bam}
        
        # Cleanup
        rm -rf $TMPIDX
        """

rule nanopore_binning_bam:
    input:
        contigs = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contigs.fasta"),
        reads = lambda wc: NANOPORE_SAMPLES[wc.sample]["fq"]
    output:
        bam = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "contigs_sorted.bam"),
        bai = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "contigs_sorted.bam.bai")
    conda:
        os.path.join(ENVS_DIR, "coverage_nanopore.yaml")
    threads: THREADS
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {output.bam})
        
        # Map reads and create sorted BAM
        minimap2 -ax map-ont -t {threads} {input.contigs} {input.reads} 2>/dev/null | \
          samtools view -bS -q 30 - | \
          samtools sort -@ {threads} -o {output.bam}
        
        # Index BAM file
        samtools index {output.bam}
        """

rule illumina_binning_depth:
    input:
        bam = os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "contigs_sorted.bam"),
        bai = os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "contigs_sorted.bam.bai")
    output:
        depth = os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "binning_depth.txt")
    conda:
        os.path.join(ENVS_DIR, "metabat2.yaml")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.depth})
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        """

rule nanopore_binning_depth:
    input:
        bam = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "contigs_sorted.bam"),
        bai = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "contigs_sorted.bam.bai")
    output:
        depth = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "binning_depth.txt")
    conda:
        os.path.join(ENVS_DIR, "metabat2.yaml")
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.depth})
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {input.bam}
        """

rule nanopore_estimate_contig_mean:
    input:
        contigs = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contigs.fasta")
    output:
        meanlen = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contig_meanlen.txt")
    conda:
        os.path.join(ENVS_DIR, "base_utils.yaml")
    shell:
        "python3 scripts/estimate_contig_lengths.py --fasta {input.contigs} --out {output.meanlen}"

rule run_metabat2_illumina:
    input:
        contigs = os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "contigs.fasta"),
        depth = os.path.join(RESULTS_DIR, "illumina", "{sample}", "qc", "binning_depth.txt")
    output:
        bins_dir = directory(os.path.join(RESULTS_DIR, "illumina", "{sample}", "binning"))
    conda:
        os.path.join(ENVS_DIR, "metabat2.yaml")
    threads: THREADS
    params:
        min_contig = lambda wildcards: config["binning"]["metabat2_min_contig"],
        min_cv = lambda wildcards: config["binning"]["metabat2_min_cv"],
        min_cv_sum = lambda wildcards: config["binning"]["metabat2_min_cv_sum"],
        max_p = lambda wildcards: config["binning"]["metabat2_max_p"],
        min_s = lambda wildcards: config["binning"]["metabat2_min_s"]
    shell:
        r"""
        mkdir -p {output.bins_dir}
        ( metabat2 \
          -i {input.contigs} \
          -a {input.depth} \
          -o {output.bins_dir}/bin \
          -m {params.min_contig} \
          -t {threads} \
          --minCV {params.min_cv} \
          --minCVSum {params.min_cv_sum} \
          --maxP {params.max_p} \
          --minS {params.min_s} \
          --verbose
        ) || true
        if ls {output.bins_dir}/bin.*.fa 2>/dev/null | grep -q .; then
          echo "MetaBAT2 bins found."
        else
          echo ">contig_unbinned" > {output.bins_dir}/bin.unbinned.fa
          echo "AAAA" >> {output.bins_dir}/bin.unbinned.fa
        fi
        """

rule run_metabat2_nanopore:
    input:
        contigs = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contigs.fasta"),
        depth = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "qc", "binning_depth.txt"),
        meanlen = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contig_meanlen.txt")
    output:
        bins_dir = directory(os.path.join(RESULTS_DIR, "nanopore", "{sample}", "binning"))
    conda:
        os.path.join(ENVS_DIR, "metabat2.yaml")
    threads: THREADS
    params:
        min_contig = lambda wildcards, input: max(1000, int(open(input.meanlen).read().strip())),
        min_cv = lambda wildcards: config["binning"]["metabat2_min_cv"],
        min_cv_sum = lambda wildcards: config["binning"]["metabat2_min_cv_sum"],
        max_p = lambda wildcards: config["binning"]["metabat2_max_p"],
        min_s = lambda wildcards: config["binning"]["metabat2_min_s"]
    shell:
        r"""
        mkdir -p {output.bins_dir}
        ( metabat2 \
          -i {input.contigs} \
          -a {input.depth} \
          -o {output.bins_dir}/bin \
          -m {params.min_contig} \
          -t {threads} \
          --minCV {params.min_cv} \
          --minCVSum {params.min_cv_sum} \
          --maxP {params.max_p} \
          --minS {params.min_s} \
          --verbose
        ) || true
        if ls {output.bins_dir}/bin.*.fa 2>/dev/null | grep -q .; then
          echo "MetaBAT2 bins found."
        else
          echo ">contig_unbinned" > {output.bins_dir}/bin.unbinned.fa
          echo "AAAA" >> {output.bins_dir}/bin.unbinned.fa
        fi
        """
