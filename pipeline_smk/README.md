# Respiratory Metagenomics Clinical Pipeline

## Pipeline Summary
This pipeline is designed for clinical metagenomic analysis of respiratory samples. It provides an end-to-end automated workflow that accepts raw sequencing data (Illumina paired-end or Nanopore single-end), performs quality control, removes host DNA, and conducts deep taxonomic profiling and antimicrobial resistance (AMR) detection.

Key features include:
- **Host Depletion**: Removes human DNA reads using Bowtie2 (Illumina) or Minimap2 (Nanopore).
- **Taxonomic Profiling**: Uses Kraken2 and Centrifuge to report precise species abundances, along with Bracken for Bayesian re-estimation.
- **AMR Detection**: Employs KMA (reads-based), AMRFinderPlus, and ABRicate (assembly-based) to robustly identify resistance genes.
- **Plasmid Reconstruction**: Uses MOB-suite to detect, type, and assess the mobility of plasmids.
- **Reporting**: Generates a comprehensive HTML and clinician-friendly PDF report detailing QC metrics, pathogen detections, and AMR profiles.

---

## Required Software Prerequisites
Ensure the following tools are installed in your environment (the pipeline uses Conda internally to manage most, but you should have Snakemake available):
- **Snakemake** (v7+ recommended)
- **Conda / Mamba**
- **Bowtie2** (for Illumina host removal)
- **Minimap2** (for Nanopore host removal)
- **KMA**
- **Kraken2 & Bracken**
- **Centrifuge**
- **AMRFinderPlus & ABRicate**
- **MOB-suite**
- **Aria2c / Curl** (for fetching external accessions)

---

## Required Databases
The pipeline requires several external databases to function. The paths to these databases must be defined in your `config/config.yaml` file.

1. **Host Genome Index**: 
   - `GRCh38` Bowtie2 index (for Illumina).
   - `GRCh38` Minimap2 `.mmi` index (for Nanopore).
2. **Kraken2 Database**: A standard or custom Kraken2 database containing archaea, bacteria, viral, and human sequences.
3. **Centrifuge Database**: E.g., `p+h+v` (prokaryotes, human, viruses).
4. **ResFinder Database**: Indexed for KMA.
5. **AMRFinderPlus Database**: Downloaded via `amrfinder_update`.
6. **MOB-suite Database**: Initialized via `mob_init`.

---

## Expected Directory Structure for Input Data
The pipeline can ingest data from URLs, SRA/ENA accessions, or local folders. 

If using a **local directory** as input, you can organize your FASTQ files freely. The ingestion script will automatically deduce the sample name based on the file naming (e.g., `_R1`/`_R2`) or use the parent folder's name.

Example acceptable structures:
```text
# Standard Illumina naming:
data/
  sampleA_R1.fastq.gz
  sampleA_R2.fastq.gz

# Custom folder naming (the script uses the folder name "Patient1"):
data/
  Patient1/
    reads_1.fq.gz
    reads_2.fq.gz
```

---

## Configuration & Results Directory
Before running, you **must** configure the pipeline using `config/config.yaml`.
Crucially, you must set the `project_root` or results directory path where the analysis outputs will be stored. Ensure you have write permissions to this path.

Example `config.yaml` snippet:
```yaml
project_root: "/path/to/my_analysis_results"
platform: "illumina" # or nanopore
```

---

## How to Run the Pipeline

The pipeline is operated via the `pipeline` wrapper script in the root directory.

### Step 1: Ingest Data
Use the `ingest` command to stage your inputs and generate a `manifests/samples.tsv` file.
```bash
./pipeline ingest
```
This will open an interactive menu allowing you to:
1. Paste URLs to FASTQ files.
2. Provide a text file with links.
3. Paste accession IDs (ERR/SRR/DRR/ERX/DRX) for automatic ENA download.
4. Point to a local directory containing `.fastq` or `.fq` files.

### Step 2: Run the Workflow
Once the manifest is generated, execute the Snakemake workflow:
```bash
./pipeline run --cores 16
```
This will automatically parse the manifest, submit jobs to Conda environments, and execute the full clinical pipeline.

---

## Summary of Outputs

Outputs are saved in your configured results directory, structured by platform and sample ID (e.g., `results/illumina/SampleA/`).

Key output directories include:
- `host_removed/`: Contains the non-host `.fastq.gz` files and a `host_stats.json` file detailing the depletion rates.
- `qc/`: Pre- and post-trimming FastQC reports and MultiQC summaries.
- `taxonomy/` & `bin_taxonomy/`: Raw Kraken2, Bracken, and Centrifuge reports.
- `assembly/`: SPAdes or Flye contigs and Quast metrics.
- `plasmids/`: MOB-suite reconstruction files (`plasmid_unbinned.fasta`, `mobtyper_results.txt`).
- `amr/`: Tabular outputs from AMRFinder, ABRicate, and KMA.
- **`reports/`**: The final **Clinical PDF and HTML reports** containing a synthesized, readable overview of the sample's metrics, detected pathogens (aggregated to species level), and combined AMR profile with coverages.
