# NGS Automated Variant Calling Pipeline

An automated Bash pipeline for germline variant detection in *E. coli* using public SRA data (`SRR1553607`).

## Prerequisites
You need [Conda](https://conda.io) or [Mamba](https://readthedocs.io) installed.

## Getting Started

1. Clone this repository:
   ```bash
   git clone https://github.com
   cd ngs-variant-calling
   ```

2. Create and activate the computational environment:
   ```bash
   conda env create -f environment.yml
   conda activate ngs-env
   ```

3. Make the script executable and run the experiment:
   ```bash
   chmod +x run_pipeline.sh
   ./run_pipeline.sh
   ```

## Pipeline Architecture
- **SRA-Tools**: Data acquisition.
- **FastQC**: Pre-alignment data assessment.
- **Trimmomatic**: Low-quality base filtering.
- **BWA-MEM**: Reference genome registration.
- **Samtools**: Alignment formatting and sorting.
- **BCFtools**: Statistical variant/mutation verification.
