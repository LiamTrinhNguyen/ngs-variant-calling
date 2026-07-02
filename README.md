## NGS Automated Variant Calling Pipeline

**An interactive Python-based Next-Generation Sequencing (NGS) pipeline for variant detection in *Escherichia coli*.**

### Input Modes

1. **Public SRA Data (Default)**: Automatically fetches paired-end reads from the NCBI Sequence Read Archive (SRA) using `fasterq-dump`.  
   *(Default sample: SRR1553607)*

2. **Custom Local FASTQ**: Supports user-provided paired-end files (`_1.fastq` and `_2.fastq`) placed in the `data/` directory.

---

### Biological & Analytical Rationale

This pipeline follows rigorous microbial NGS best practices:

- **Integrated QC & Trimming**: Uses **fastp** for high-performance adapter trimming, quality filtering, and base correction in a single efficient pass.
- **Statistical Precision**: Employs **BCFtools** for variant calling, combined with an automated **Polars** filtering layer that removes low-confidence calls (DP < 15, QUAL ≤ 30).
- **Functional Annotation**: Integrates **SnpEff** to predict the biological impact of variants (e.g., missense, synonymous, frameshift, stop-gain).
- **Resource Management**: Automatic cleanup of large intermediate files (SAM, raw FASTQs) to maintain a minimal storage footprint.

---

### Prerequisites

- **Environment**: Linux (recommended: GitHub Codespaces or Ubuntu)
- **Package Manager**: Conda or Mamba
- **Core Dependencies**:
  - `fastp`, `bwa`, `samtools`, `bcftools`, `snpEff` (v5.1+)
  - OpenJDK 11
  - Python packages: `polars`, `tqdm`
- **Storage**: Approximately **2–4 GB** of free disk space (optimized pipeline)

---

### Pipeline Steps Overview

| Step | Tool                  | Purpose                                      | Biological Reasoning                                      |
|------|-----------------------|----------------------------------------------|-----------------------------------------------------------|
| 1    | Setup                 | Create organized directories                 | Maintains reproducibility and prevents data mixing        |
| 2    | Reference             | Download *E. coli* K-12 reference genome     | Provides a trusted baseline for accurate comparison       |
| 3    | Indexing              | Build BWA and samtools indices               | Enables fast and accurate read mapping                    |
| 4    | Data Input            | Download from SRA or use local FASTQ         | Flexibility for public or private datasets                |
| 5    | fastp                 | QC, trimming, and base correction            | Removes noise/adapters; generates QC metrics              |
| 6    | BWA-MEM + Samtools    | Align reads and prepare BAM files            | Accurate mapping to reference genome                      |
| 7    | BCFtools              | Variant calling                              | Statistical identification of true mutations              |
| 8    | SnpEff                | Functional annotation                        | Predicts biological impact of variants                    |
| 9    | Polars                | Data filtering & reporting                   | Extracts high-confidence signals from VCF noise           |

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
   - final_variants.vcf – Raw detected mutations.
   - annotated_variants.vcf – Functional impact predictions (via SnpEff).
   - fastp.html/json – Comprehensive quality and adapter content reports.
   - alignment_stats.txt – Alignment mapping rates for data validation.
   - pipeline_*.log – Full audit trail including the final filtered mutation table.

## Features
   - Automated Validation: Step 8 logs the alignment rate and warns if coverage is insufficient.
   - Data Decomposition: Uses polars to extract Depth (DP), Allele Count (AC), and Effect into a clean, filterable table.
   - Production-Ready: Includes error handling, auto-cleanup of intermediate SAM/FASTQ files, and Java memory allocation (-Xmx4g) for annotation.
   - Interactive Dashboard: Single-step initial confirmation followed by fully automated execution.

## Tools Used
   - Data Retrieval: SRA-Tools (fasterq-dump)
   - QC & Trimming: fastp
   - Contamination Check: Kraken2
   - Alignment: BWA-MEM + Samtools
   - Variant Calling: BCFtools
   - Annotation: SnpEff
   - Data Processing: Polars + tqdm

## Author
Liam TrinhNguyen
   Data Scientist, Wisconsin State Laboratory of Hygiene (WSLH)
   Bioinformatics and Computational Biology, M.S.
