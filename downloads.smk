import os
import re
import subprocess

def _aria2c_common_args():
    a = config["download"]["aria2c"]
    args = [
        "aria2c",
        "-x", str(a.get("connections", 16)),
        "-s", str(a.get("split", 16)),
        "-j", str(a.get("parallel", 4)),
        "--retry-wait", str(a.get("retry_wait", 5)),
        "--max-tries", str(a.get("max_tries", 5)),
    ]
    if a.get("check_integrity", True):
        args.append("--check-integrity=true")
    return args

def _ena_fastq_base_url(acc):
    prefix = acc[:6]
    subdir = ""
    if len(acc) > 9:
        try:
            subdir = "/" + f"{int(acc[9:]):03d}"
        except:
            subdir = ""
    return f"ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{prefix}{subdir}/{acc}"

# ------------- ILLUMINA ----------------
rule download_illumina_sra_r1:
    output:
        fq=lambda wc: os.path.join(PROJECT_ROOT, "illumina", wc.acc, f"{wc.acc}_1.fastq.gz"),
        done=lambda wc: os.path.join(PROJECT_ROOT, "illumina", wc.acc, ".download.done"),
        valid=lambda wc: os.path.join(PROJECT_ROOT, "illumina", wc.acc, ".r1.validated")
    params:
        url=lambda wc: f"{_ena_fastq_base_url(wc.acc)}/{wc.acc}_1.fastq.gz"
    conda:
        "envs/aria2c.yaml"
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.fq})
        {" ".join(_aria2c_common_args())} -d $(dirname {output.fq}) "{params.url}"
        touch {output.done}
        python3 scripts/validate_fastq.py --fastq {output.fq}
        touch {output.valid}
        """

rule download_illumina_sra_r2:
    output:
        fq=lambda wc: os.path.join(PROJECT_ROOT, "illumina", wc.acc, f"{wc.acc}_2.fastq.gz"),
        done=lambda wc: os.path.join(PROJECT_ROOT, "illumina", wc.acc, ".download.done"),
        valid=lambda wc: os.path.join(PROJECT_ROOT, "illumina", wc.acc, ".r2.validated")
    params:
        url=lambda wc: f"{_ena_fastq_base_url(wc.acc)}/{wc.acc}_2.fastq.gz"
    conda:
        "envs/aria2c.yaml"
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.fq})
        {" ".join(_aria2c_common_args())} -d $(dirname {output.fq}) "{params.url}"
        touch {output.done}
        python3 scripts/validate_fastq.py --fastq {output.fq}
        touch {output.valid}
        """

# ------------- NANOPORE ----------------
rule download_nanopore_sra:
    output:
        fq=lambda wc: os.path.join(PROJECT_ROOT, "nanopore", wc.acc, f"{wc.acc}_1.fastq.gz"),
        done=lambda wc: os.path.join(PROJECT_ROOT, "nanopore", wc.acc, ".download.done"),
        valid=lambda wc: os.path.join(PROJECT_ROOT, "nanopore", wc.acc, ".r1.validated")
    params:
        url1=lambda wc: f"{_ena_fastq_base_url(wc.acc)}/{wc.acc}_1.fastq.gz",
        url2=lambda wc: f"{_ena_fastq_base_url(wc.acc)}/{wc.acc}.fastq.gz"
    conda:
        "envs/aria2c.yaml"
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.fq})
        # Try _1.fastq.gz first, fallback to .fastq.gz
        {" ".join(_aria2c_common_args())} -d $(dirname {output.fq}) "{params.url1}" || \
        {" ".join(_aria2c_common_args())} -d $(dirname {output.fq}) "{params.url2}"
        touch {output.done}
        python3 scripts/validate_fastq.py --fastq {output.fq}
        touch {output.valid}
        """

# ------------- GENERIC URL -------------
rule download_url_generic:
    output:
        fq=os.path.join(PROJECT_ROOT, "{platform}", "url_downloads", "{fname}"),
        done=os.path.join(PROJECT_ROOT, "{platform}", "url_downloads", ".download.done"),
        valid=os.path.join(PROJECT_ROOT, "{platform}", "url_downloads", "{fname}.validated")
    params:
        url=lambda wc: wc.url
    conda:
        "envs/aria2c.yaml"
    threads: 1
    shell:
        r"""
        mkdir -p $(dirname {output.fq})
        {" ".join(_aria2c_common_args())} -d $(dirname {output.fq}) "{params.url}"
        touch {output.done}
        python3 scripts/validate_fastq.py --fastq {output.fq}
        touch {output.valid}
        """
