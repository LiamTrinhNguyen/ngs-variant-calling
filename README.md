# NGS Automated Variant Calling Pipeline

An interactive Python-based Next-Generation Sequencing (NGS) pipeline for **variant detection** in *Escherichia coli*.

The pipeline supports **two input modes**:
1. **Public SRA data** (default) – automatically downloads from NCBI SRA.
2. **Custom local FASTQ files** – you can provide your own paired-end FASTQ files.

---

## Biological & Analytical Rationale

This pipeline follows standard microbial NGS best practices to transform raw sequencing data into reliable variant calls. It emphasizes quality control, accurate alignment, and statistical variant detection to minimize false positives while maintaining biological relevance.

---

## Prerequisites

- **Conda** or **Mamba**
- Internet access (required only for public SRA mode)
- ~2–10 GB of disk space

---


## Getting Started

1. **Clone the repository**
   ```bash
   git clone https://github.com/yourusername/ngs-variant-calling.git
   cd ngs-variant-calling
   ```
2. **Create and activate the environment**
   ```Bash
   conda env create -f environment.yml
   conda activate ngs-env
   ```
   
3. **Run the pipeline**
   ```Bash
   python3 run_interactive_pipeline.py
   ```






## Pipeline Workflow

| Step | Tool                    | Purpose                                           | Biological Reasoning                                              |
|------|-------------------------|---------------------------------------------------|-------------------------------------------------------------------|
| 1    | Setup                   | Create organized directories                      | Maintains reproducibility and clear project structure             |
| 2    | Reference               | Download *E. coli* K-12 reference genome          | Provides trusted baseline for variant comparison                  |
| 3    | Indexing                | Build BWA and samtools indices                    | Enables fast and accurate read mapping                            |
| 4    | Data Input              | Retrieve SRA data or use local FASTQ              | Flexibility for public or private sequencing datasets             |
| 5    | **fastp**               | Combined QC, trimming, and base correction        | Removes noise/adapters and generates comprehensive QC report      |
| 6    | BWA-MEM + Samtools      | Align reads and prepare sorted/indexed BAM        | Accurate mapping of reads to reference genome                     |
| 7    | BCFtools                | Variant calling                                   | Statistical detection of mutations and indels                     |
| 8    | SnpEff                  | Functional annotation of variants                 | Predicts biological consequences (e.g., missense, stop-gain)     |
| 9    | Polars + Reporting      | Filtering & final mutation report                 | Extracts high-quality variants and presents clean results        |





## Input Options

   1. **Option1: Public SRA Data (Default)**
   The pipeline automatically downloads paired-end FASTQ files using fasterq-dump.
   Default sample: SRR1553607
   2. **Option2: Custom Local FASTQ Files**
   
   Place your paired-end FASTQ files (.fastq or .fastq.gz) in the data/ folder (or any accessible location).
   When prompted, select Option 2 and enter the base name of your files (without _1 / _2 and extensions).
   
   Example:
   Files: my_sample_1.fastq.gz and my_sample_2.fastq.gz -> Enter: my_sample

## Outputs
   Each run creates a timestamped folder under results/run_YYYYMMDD_HHMMSS/:
   final_variants.vcf – Raw detected mutations.
   annotated_variants.vcf – Functional impact predictions (via SnpEff).
   fastp.html/json – Comprehensive quality and adapter content reports.
   alignment_stats.txt – Alignment mapping rates for data validation.
   pipeline_*.log – Full audit trail including the final filtered mutation table.

## Features
   Automated Validation: Step 8 logs the alignment rate and warns if coverage is insufficient.
   Data Decomposition: Uses polars to extract Depth (DP), Allele Count (AC), and Effect into a clean, filterable table.
   Production-Ready: Includes error handling, auto-cleanup of intermediate SAM/FASTQ files, and Java memory allocation (-Xmx4g) for annotation.
   Interactive Dashboard: Single-step initial confirmation followed by fully automated execution.

## Tools Used
   Data Retrieval: SRA-Tools (fasterq-dump)
   QC & Trimming: fastp
   Contamination Check: Kraken2
   Alignment: BWA-MEM + Samtools
   Variant Calling: BCFtools
   Annotation: SnpEff
   Data Processing: Polars + tqdm

## Author
Liam TrinhNguyen
   Data Scientist, Wisconsin State Laboratory of Hygiene (WSLH)
   Bioinformatics and Computational Biology, M.S.
