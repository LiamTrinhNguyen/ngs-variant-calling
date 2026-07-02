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

| Step | Tool                    | Purpose                                      | Biological Reasoning |
|------|-------------------------|----------------------------------------------|----------------------|
| 1    | Setup                   | Create organized directories                 | Maintains reproducibility |
| 2    | Reference               | Download *E. coli* reference genome          | Provides trusted baseline for comparison |
| 3    | Indexing                | Build alignment indices                      | Enables fast and accurate read mapping |
| 4    | Data Input              | Download from SRA **or** use local FASTQ files | Flexibility for public or private datasets |
| 5    | FastQC                  | Quality assessment                           | Detects biases and quality issues |
| 6    | Trimmomatic             | Trim adapters and low-quality bases          | Removes noise that could cause false variants |
| 7    | BWA-MEM + Samtools      | Align reads and prepare BAM files            | Accurate mapping to reference genome |
| 8    | BCFtools                | Variant calling                              | Statistical identification of true mutations |





## Input Options

1. **Option1: Public SRA Data (Default)**
The pipeline automatically downloads paired-end FASTQ files using fasterq-dump.
Default sample: SRR1553607
2. **Option2: Custom Local FASTQ Files**

Place your paired-end FASTQ files (.fastq or .fastq.gz) in the data/ folder (or any accessible location).
When prompted, select Option 2 and enter the base name of your files (without _1 / _2 and extensions).

Example:
Files: my_sample_1.fastq.gz and my_sample_2.fastq.gz → Enter: my_sample

## Outputs
   Each run creates a timestamped folder under results/run_YYYYMMDD_HHMMSS/ to prevent overwriting previous results:
   
   final_variants.vcf – List of detected mutations
   FastQC quality reports
   Sorted & indexed BAM file
   Detailed log file (including the full Polars mutation table)


## Features

   Only Step 1 requires manual confirmation; remaining steps run automatically
   Real-time progress bars (tqdm)
   Timestamped run directories (safe for repeated runs)
   Comprehensive logging with complete mutation table
   Support for both public SRA and custom local FASTQ inputs


## Tools Used

   Data Retrieval: SRA-Tools (fasterq-dump)
   Quality Control: FastQC
   Trimming: Trimmomatic
   Alignment: BWA-MEM + Samtools
   Variant Calling: BCFtools
   Reporting: Polars + tqdm


## Author

   Liam TrinhNguyen
   APHL-CDC Public Health Laboratory Fellow / Data Scientist
   Wisconsin State Laboratory of Hygiene
