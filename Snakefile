import os
import csv
import json
from pathlib import Path

configfile: "config/config.yaml"

PROJECT_ROOT = os.path.realpath(os.path.join(workflow.basedir, config.get("project_root", "..")))
RESULTS_DIR  = os.path.normpath(config.get("results_dir", "results"))
ENVS_DIR = os.path.join(workflow.basedir, "envs")

PROFILE  = (config.get("profile", "bioinformatician") or "bioinformatician").strip()
PLATFORM = (config.get("platform", "illumina") or "illumina").strip().lower()
THREADS  = int(config.get("threads", 8))

MANIFEST = os.path.join(PROJECT_ROOT, "manifests", "samples.tsv")
DB_JSON  = os.path.join(RESULTS_DIR, "detected_databases.json")


def read_manifest():
    illum = {}
    nano = {}
    if not os.path.exists(MANIFEST):
        return illum, nano
    with open(MANIFEST, "r") as fh:
        raw_lines = [line.rstrip("\n") for line in fh if line.strip() and not line.lstrip().startswith("#")]
    if not raw_lines:
        return illum, nano
    header = raw_lines[0]
    if "\t" in header:
        cols = [c.strip() for c in header.split("\t")]
        use_tab = True
    else:
        cols = [c.strip() for c in header.split()]
        use_tab = False
    idx = {c.lower(): i for i, c in enumerate(cols)}
    def split_line(line):
        if use_tab and "\t" in line:
            return line.split("\t")
        return line.split()
    for line in raw_lines[1:]:
        parts = split_line(line)
        if not parts:
            continue
        platform = ""
        sample = ""
        if "platform" in idx and idx["platform"] < len(parts):
            platform = parts[idx["platform"]].strip().lower()
        elif len(parts) > 0:
            platform = parts[0].strip().lower()
        if "sample" in idx and idx["sample"] < len(parts):
            sample = parts[idx["sample"]].strip()
        elif len(parts) > 1:
            sample = parts[1].strip()
        if not sample:
            continue
        if platform.startswith("illum"):
            if "r1" in idx and idx["r1"] < len(parts):
                r1 = parts[idx["r1"]].strip()
            elif len(parts) > 2:
                r1 = parts[2].strip()
            else:
                r1 = ""
            if "r2" in idx and idx["r2"] < len(parts):
                r2 = parts[idx["r2"]].strip()
            elif len(parts) > 3:
                r2 = parts[3].strip()
            else:
                r2 = ""
            if r1:
                illum[sample] = {"r1": r1, "r2": r2}
        else:
            if "fq" in idx and idx["fq"] < len(parts):
                fq = parts[idx["fq"]].strip()
            elif len(parts) > 2:
                fq = parts[2].strip()
            else:
                fq = ""
            if fq:
                nano[sample] = {"fq": fq}
    return illum, nano

ILLUMINA_SAMPLES, NANOPORE_SAMPLES = read_manifest()
SAMPLES_ILLUMINA = sorted(ILLUMINA_SAMPLES.keys())
SAMPLES_NANOPORE = sorted(NANOPORE_SAMPLES.keys())

if PLATFORM not in ("illumina", "nanopore"):
    raise ValueError(f"Invalid config platform: {PLATFORM!r}. Allowed: illumina | nanopore")

if not os.path.exists(MANIFEST):
    raise ValueError(
        f"Manifest not found: {MANIFEST}\n"
        f"Run: ./pipeline ingest (or: ./pipeline run --interactive-manifest)"
    )

if not SAMPLES_ILLUMINA and not SAMPLES_NANOPORE:
    raise ValueError(
        f"Manifest has no usable sample rows: {MANIFEST}\n"
        f"Expected columns: platform, sample, r1, r2 (tab-separated) or platform sample fq for nanopore.\n"
        f"For illumina rows, r1 must be present; r2 is expected for paired-end."
    )

if PLATFORM == "illumina" and not SAMPLES_ILLUMINA:
    raise ValueError(
        f"Config platform is 'illumina' but no illumina samples were found in manifest: {MANIFEST}\n"
        f"Either ingest illumina paired reads or set platform: 'nanopore'."
    )

if PLATFORM == "nanopore" and not SAMPLES_NANOPORE:
    raise ValueError(
        f"Config platform is 'nanopore' but no nanopore samples were found in manifest: {MANIFEST}\n"
        f"Either ingest nanopore reads or set platform: 'illumina'."
    )

def find_up_rglob(start_dir, pattern, max_parent_levels=6):
    cur = Path(start_dir).resolve()
    for _ in range(max_parent_levels):
        hits = list(cur.rglob(pattern))
        hits = [h for h in hits if h.is_file()]
        if hits:
            return str(hits[0])
        if cur.parent == cur:
            break
        cur = cur.parent
    return None

def detect_kraken_db():
    h = find_up_rglob(PROJECT_ROOT, "hash.k2d")
    return os.path.dirname(h) if h else None

def detect_centrifuge_prefix():
    cf = find_up_rglob(PROJECT_ROOT, "*.1.cf")
    return str(cf).replace(".1.cf", "") if cf else None

def detect_host_mmi():
    return find_up_rglob(PROJECT_ROOT, "GRCh38*.mmi")

def detect_host_bt2_prefix():
    bt2 = find_up_rglob(PROJECT_ROOT, "*.1.bt2")
    return str(bt2).replace(".1.bt2", "") if bt2 else None

def detect_checkm2_db():
    preferred = Path(PROJECT_ROOT) / "checkm2_db" / "CheckM2_database" / "CheckM2_database" / "uniref100.KO.1.dmnd"
    if preferred.is_file():
        return str(preferred)
    return find_up_rglob(PROJECT_ROOT, "uniref100.KO.1.dmnd")

def detect_resfinder_kma_prefix():
    for base in [Path(PROJECT_ROOT), Path(PROJECT_ROOT).parent]:
        for p in base.rglob("resfinder_db"):
            seqb = list(Path(p).rglob("resfinder_db_kma.seq.b"))
            if seqb:
                return str(seqb[0]).replace(".seq.b", "")
    return None

def illumina_reports():
    if PROFILE == "clinician":
        return [os.path.join(RESULTS_DIR, "illumina", s, "report", "report_clinician.json") for s in SAMPLES_ILLUMINA]
    else:
        return [os.path.join(RESULTS_DIR, "illumina", s, "report", "report_bioinformatician.json") for s in SAMPLES_ILLUMINA]

def illumina_reports_txt():
    if PROFILE == "clinician":
        return [os.path.join(RESULTS_DIR, "illumina", s, "report", "report_clinician.txt") for s in SAMPLES_ILLUMINA]
    else:
        return [os.path.join(RESULTS_DIR, "illumina", s, "report", "report_bioinformatician.txt") for s in SAMPLES_ILLUMINA]

def nanopore_reports():
    if PROFILE == "clinician":
        return [os.path.join(RESULTS_DIR, "nanopore", s, "report", "report_clinician.json") for s in SAMPLES_NANOPORE]
    else:
        return [os.path.join(RESULTS_DIR, "nanopore", s, "report", "report_bioinformatician.json") for s in SAMPLES_NANOPORE]

def nanopore_reports_txt():
    if PROFILE == "clinician":
        return [os.path.join(RESULTS_DIR, "nanopore", s, "report", "report_clinician.txt") for s in SAMPLES_NANOPORE]
    else:
        return [os.path.join(RESULTS_DIR, "nanopore", s, "report", "report_bioinformatician.txt") for s in SAMPLES_NANOPORE]

def illumina_centrifuge_targets():
    outs = []
    for s in SAMPLES_ILLUMINA:
        outs += [
            os.path.join(RESULTS_DIR, "illumina", s, "taxonomy", "centrifuge_output.txt"),
            os.path.join(RESULTS_DIR, "illumina", s, "taxonomy", "centrifuge_report.tsv"),
            os.path.join(RESULTS_DIR, "illumina", s, "taxonomy", "centrifuge_kreport.txt"),
        ]
    return outs

def nanopore_centrifuge_targets():
    outs = []
    for s in SAMPLES_NANOPORE:
        outs += [
            os.path.join(RESULTS_DIR, "nanopore", s, "taxonomy", "centrifuge_output.txt"),
            os.path.join(RESULTS_DIR, "nanopore", s, "taxonomy", "centrifuge_report.tsv"),
            os.path.join(RESULTS_DIR, "nanopore", s, "taxonomy", "centrifuge_kreport.txt"),
        ]
    return outs

def illumina_krona_targets():
    return [os.path.join(RESULTS_DIR, "illumina", s, "taxonomy", "krona_kraken2.html") for s in SAMPLES_ILLUMINA] + \
           [os.path.join(RESULTS_DIR, "illumina", s, "taxonomy", "krona_centrifuge.html") for s in SAMPLES_ILLUMINA]

def nanopore_krona_targets():
    return [os.path.join(RESULTS_DIR, "nanopore", s, "taxonomy", "krona_kraken2.html") for s in SAMPLES_NANOPORE] + \
           [os.path.join(RESULTS_DIR, "nanopore", s, "taxonomy", "krona_centrifuge.html") for s in SAMPLES_NANOPORE]

def clinician_extra_amr_targets():
    outs = []
    for s in SAMPLES_ILLUMINA:
        outs += [
            os.path.join(RESULTS_DIR, "illumina", s, "amr", "abricate_card.tsv"),
            os.path.join(RESULTS_DIR, "illumina", s, "amr", "amrfinder_results.tsv"),
            os.path.join(RESULTS_DIR, "illumina", s, "amr", "kma_resfinder", "kma.res"),
            os.path.join(RESULTS_DIR, "illumina", s, "amr", "kma_resfinder", "kma.tsv"),
        ]
    for s in SAMPLES_NANOPORE:
        outs += [
            os.path.join(RESULTS_DIR, "nanopore", s, "amr", "abricate_card.tsv"),
            os.path.join(RESULTS_DIR, "nanopore", s, "amr", "amrfinder_results.tsv"),
        ]
    return outs

rule write_detected_dbs:
    input:
        MANIFEST
    output:
        DB_JSON
    run:
        os.makedirs(RESULTS_DIR, exist_ok=True)
        data = {
            "project_root": PROJECT_ROOT,
            "manifest": MANIFEST,
            "kraken_db": detect_kraken_db(),
            "centrifuge_prefix": detect_centrifuge_prefix(),
            "host_mmi": detect_host_mmi(),
            "host_bt2_prefix": detect_host_bt2_prefix(),
            "checkm2_db": detect_checkm2_db(),
            "resfinder_kma_prefix": detect_resfinder_kma_prefix(),
            "samples_detected": {"illumina": SAMPLES_ILLUMINA, "nanopore": SAMPLES_NANOPORE},
            "profile": PROFILE,
            "platform": PLATFORM
        }
        with open(output[0], "w") as f:
            json.dump(data, f, indent=2)

rule all:
    input:
        DB_JSON,
        illumina_centrifuge_targets() if PLATFORM == "illumina" else [],
        nanopore_centrifuge_targets() if PLATFORM == "nanopore" else [],
        illumina_krona_targets() if PLATFORM == "illumina" else [],
        nanopore_krona_targets() if PLATFORM == "nanopore" else [],
        illumina_reports() if PLATFORM == "illumina" else [],
        illumina_reports_txt() if PLATFORM == "illumina" else [],
        nanopore_reports() if PLATFORM == "nanopore" else [],
        nanopore_reports_txt() if PLATFORM == "nanopore" else [],
        clinician_extra_amr_targets() if PROFILE == "clinician" else []

include: "rules/common.smk"
include: "rules/qc_illumina.smk"
include: "rules/host_removal_illumina.smk"
include: "rules/taxonomy_illumina.smk"
include: "rules/qc_nanopore.smk"
include: "rules/host_removal_nanopore.smk"
include: "rules/taxonomy_nanopore.smk"
include: "rules/illumina.smk"
include: "rules/nanopore.smk"
include: "rules/assembly.smk"
include: "rules/quast.smk"
include: "rules/coverage.smk"
include: "rules/amr.smk"
include: "rules/assembly_nanopore.smk"
include: "rules/quast_nanopore.smk"
include: "rules/coverage_nanopore.smk"
include: "rules/amr_nanopore.smk"
include: "rules/binning.smk"
include: "rules/checkm2.smk"
include: "rules/prokka.smk"
include: "rules/bin_amr_illumina.smk"
include: "rules/bin_amr_nanopore.smk"
