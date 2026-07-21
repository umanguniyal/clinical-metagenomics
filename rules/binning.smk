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
          --verbose \
          --unbinned
        ) || true
        if ls {output.bins_dir}/bin.*.fa 2>/dev/null | grep -q .; then
          echo "MetaBAT2 bins found."
        fi
        if [ ! -f {output.bins_dir}/bin.unbinned.fa ]; then
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
          --verbose \
          --unbinned
        ) || true
        if ls {output.bins_dir}/bin.*.fa 2>/dev/null | grep -q .; then
          echo "MetaBAT2 bins found."
        fi
        if [ ! -f {output.bins_dir}/bin.unbinned.fa ]; then
          echo ">contig_unbinned" > {output.bins_dir}/bin.unbinned.fa
          echo "AAAA" >> {output.bins_dir}/bin.unbinned.fa
        fi
        """

rule run_mob_recon_illumina:
    input:
        assembly = os.path.join(RESULTS_DIR, "illumina", "{sample}", "assembly", "contigs.fasta")
    output:
        outdir = directory(os.path.join(RESULTS_DIR, "illumina", "{sample}", "plasmids")),
        plasmid_fasta = os.path.join(RESULTS_DIR, "illumina", "{sample}", "plasmids", "plasmid_unbinned.fasta"),
        mobtyper = os.path.join(RESULTS_DIR, "illumina", "{sample}", "plasmids", "mobtyper_results.txt"),
        contig_report = os.path.join(RESULTS_DIR, "illumina", "{sample}", "plasmids", "contig_report.txt")
    conda:
        os.path.join(ENVS_DIR, "mob_suite.yaml")
    threads: THREADS
    shell:
        r"""
        rm -rf {output.outdir}
        mob_recon --infile {input.assembly} --outdir {output.outdir} --num_threads {threads} --force || true
        mkdir -p {output.outdir}
        touch {output.outdir}/plasmid_unbinned.fasta
        touch {output.plasmid_fasta}
        touch {output.mobtyper}
        touch {output.contig_report}
        """

rule run_mob_recon_nanopore:
    input:
        assembly = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "assembly", "contigs.fasta")
    output:
        outdir = directory(os.path.join(RESULTS_DIR, "nanopore", "{sample}", "plasmids")),
        plasmid_fasta = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "plasmids", "plasmid_unbinned.fasta"),
        mobtyper = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "plasmids", "mobtyper_results.txt"),
        contig_report = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "plasmids", "contig_report.txt")
    conda:
        os.path.join(ENVS_DIR, "mob_suite.yaml")
    threads: THREADS
    shell:
        r"""
        rm -rf {output.outdir}
        mob_recon --infile {input.assembly} --outdir {output.outdir} --num_threads {threads} --force || true
        mkdir -p {output.outdir}
        touch {output.outdir}/plasmid_unbinned.fasta
        touch {output.plasmid_fasta}
        touch {output.mobtyper}
        touch {output.contig_report}
        """

rule illumina_classify_bins:
    input:
        bins_dir = os.path.join(RESULTS_DIR, "illumina", "{sample}", "binning"),
        dbjson = DB_JSON
    output:
        tax_dir = directory(os.path.join(RESULTS_DIR, "illumina", "{sample}", "bin_taxonomy"))
    conda:
        os.path.join(ENVS_DIR, "centrifuge.yaml")
    threads: THREADS
    shell:
        r"""
        mkdir -p {output.tax_dir}
        CF_DB=$(python3 -c "import json; print(json.load(open('{input.dbjson}'))['centrifuge_prefix'])")
        for bin_fa in {input.bins_dir}/bin.*.fa; do
            [ -f "$bin_fa" ] || continue
            bin_name=$(basename "$bin_fa" .fa)
            if [ "$bin_name" != "bin.unbinned" ]; then
                centrifuge -f -x $CF_DB -U $bin_fa -S {output.tax_dir}/$bin_name.cf.out \
                    --report-file {output.tax_dir}/$bin_name.cf.report.tsv -p {threads} || true
                # Extract best hit from report
                if [ -s {output.tax_dir}/$bin_name.cf.report.tsv ]; then
                    # Skip header and get top species
                    awk -F'\t' '$3=="species" || $3=="leaf" {{print $0}}' {output.tax_dir}/$bin_name.cf.report.tsv | sort -t$'\t' -k6,6nr | head -n 1 > {output.tax_dir}/$bin_name.best_hit.txt || true
                fi
            fi
        done
        """

rule nanopore_classify_bins:
    input:
        bins_dir = os.path.join(RESULTS_DIR, "nanopore", "{sample}", "binning"),
        dbjson = DB_JSON
    output:
        tax_dir = directory(os.path.join(RESULTS_DIR, "nanopore", "{sample}", "bin_taxonomy"))
    conda:
        os.path.join(ENVS_DIR, "centrifuge.yaml")
    threads: THREADS
    shell:
        r"""
        mkdir -p {output.tax_dir}
        CF_DB=$(python3 -c "import json; print(json.load(open('{input.dbjson}'))['centrifuge_prefix'])")
        for bin_fa in {input.bins_dir}/bin.*.fa; do
            [ -f "$bin_fa" ] || continue
            bin_name=$(basename "$bin_fa" .fa)
            if [ "$bin_name" != "bin.unbinned" ]; then
                centrifuge -f -x $CF_DB -U $bin_fa -S {output.tax_dir}/$bin_name.cf.out \
                    --report-file {output.tax_dir}/$bin_name.cf.report.tsv -p {threads} || true
                if [ -s {output.tax_dir}/$bin_name.cf.report.tsv ]; then
                    awk -F'\t' '$3=="species" || $3=="leaf" {{print $0}}' {output.tax_dir}/$bin_name.cf.report.tsv | sort -t$'\t' -k6,6nr | head -n 1 > {output.tax_dir}/$bin_name.best_hit.txt || true
                fi
            fi
        done
        """
