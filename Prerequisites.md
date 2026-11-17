# Prerequisites

Everything needed to run `ngs-variant-calling/run_interactive_pipeline.py`.

## 1. Command-line tools (via Conda/Mamba, bioconda channel)

| Tool | Package | Used for |
|------|---------|----------|
| `fastp` | `fastp` | QC, adapter trimming, base correction |
| `bwa` | `bwa` | Reference indexing + BWA-MEM alignment |
| `samtools` | `samtools` | faidx, dict, sort, index |
| `bcftools` | `bcftools` | mpileup, call, norm, filter, index, query |
| `bgzip`, `tabix` | `htslib` | Compressing/indexing annotated VCF |
| `snpEff` | `snpeff` | Functional annotation |
| `java` | `openjdk` | Runtime required by SnpEff |

```bash
# Install Miniforge (skip if conda/mamba already available)
#   Miniforge is a minimal installer for the conda package manager.
#   What it is for:
#     - Provides "conda" and the faster "mamba" to install software and manage isolated environments
#     - Uses conda-forge as its default channel (free, community-maintained, no Anaconda licence limits)
#     - Most bioinformatics tools (bwa, samtools, bcftools, fastp, snpEff) ship as conda packages
#       on bioconda, and conda installs them with all their C/Java dependencies; apt or pip can't
#     - Each project gets its own environment (here "ngs-env"), so tool versions don't clash
#       with system packages or with other projects

# Download the Linux x86_64 installer script
wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh

# Run the installer:
#   -b : batch mode (no interactive questions, licence accepted automatically)
#   -p : install prefix, here ~/miniforge3
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"

# Load conda/mamba into the current shell session (activates the "base" environment)
source "$HOME/miniforge3/bin/activate"

# Optional: make conda available automatically in every new terminal
# conda init bash

# Create the environment with all tools + Python packages
#   -n ngs-env                 : environment name
#   -c conda-forge -c bioconda : channels (bioconda hosts the bioinformatics tools)
#   -y                         : skip the confirmation prompt
mamba create -n ngs-env -c conda-forge -c bioconda -y \
    python=3.11 \
    fastp bwa samtools bcftools htslib snpeff openjdk \
    polars tqdm requests biopython \
    snakemake-minimal nextflow
# A trailing "\" can't be followed by a comment, so the packages are explained here:
#   python=3.11 : interpreter that runs run_interactive_pipeline.py
#   fastp       : step 5, read QC, adapter trimming and base correction (HTML/JSON reports)
#   bwa         : step 3 "bwa index" builds the reference index; step 6 "bwa mem" aligns reads
#   samtools    : faidx/dict index the reference; sort/index turn the SAM into a sorted BAM
#   bcftools    : step 7, "mpileup -q 20 -Q 20 | call --ploidy 1 | norm" calls SNPs/indels as haploid
#                 from good-quality reads only, and writes each variant one standard way;
#                 step 9, "filter" flags calls (LowQual/LowDepth/LowAF) and "query" pulls the columns
#   htslib      : provides bgzip and tabix, which compress and index the annotated VCF
#   snpeff      : step 8, predicts variant effects (missense, frameshift, stop-gain, ...)
#   openjdk     : Java runtime that SnpEff needs
#   polars      : step 9, fast DataFrame library that computes allele frequency (AF) and writes
#                 the PASS and all-calls TSV reports
#   tqdm        : progress bars while FASTQ files download
#   requests    : HTTP client for ENA metadata and FASTQ downloads (SRR/ERR/DRR accessions)
#   biopython   : Bio.Entrez/SeqIO, which download reference FASTA and metadata from NCBI
#   snakemake-minimal : workflow manager for the multi-sample route (workflow/Snakefile).
#                 It runs the same steps for every sample in config/config.yaml, in parallel,
#                 and skips steps whose outputs already exist. ("-minimal" leaves out
#                 cloud/cluster plugins that a single machine doesn't need.)
#                 Only needed for Snakemake; the interactive script works without it.
#   nextflow    : the second workflow manager (workflow/main.nf). Same steps, same
#                 config/config.yaml; also runs on SLURM/PBS clusters and in the cloud.
#                 Uses the openjdk already listed above. Only needed for the Nextflow route.

# Already created ngs-env before Snakemake/Nextflow were added? Add them to the existing env:
#   mamba install -n ngs-env -c conda-forge -c bioconda -y snakemake-minimal nextflow

## 1b. R, in its own environment (step 9b · statistics)
# WHY a separate env: r-base needs a newer GSL library than bcftools 1.25. Installing
#     r-base INTO ngs-env makes mamba downgrade bcftools to 1.5, which then fails with
#     "libcrypto.so.1.0.0: cannot open shared object file". Two envs avoid the clash.
# Only base R is used (no CRAN packages), so r-base is all you need.
mamba create -n ngs-r -c conda-forge -y r-base
# You do NOT activate ngs-r: the pipeline finds ~/miniforge3/envs/ngs-r/bin/Rscript itself.
# To use another R, set:  export NGS_RSCRIPT=/path/to/Rscript
# Without R the interactive script still finishes (only warns); Snakemake and Nextflow stop.

## 1a.  Activate the environment (run this in every new terminal before running the pipeline)
conda activate ngs-env

# Troubleshooting:
#   "snpEff: command not found"          -> ngs-env is not active; activate it first
#   "critical libmamba Shell not initialized" when running "mamba activate"
#                                         -> mamba 2.x needs a shell hook; skip it and use:
#   source "$HOME/miniforge3/bin/activate" ngs-env
#   (or, if you want to use mamba:  eval "$(mamba shell hook --shell bash)" && mamba activate ngs-env)
```

## 2. Python packages

| Package | Required? | Used for |
|---------|-----------|----------|
| `polars` | Required (imported at startup) | Reads the flagged VCF, computes AF, writes the TSV reports |
| `tqdm` | Auto-installed if missing | Download progress bars |
| `requests` | Required for SRA/ENA accessions | FASTQ + metadata download from ENA |
| `biopython` | Required for NCBI accessions/reference download | Entrez FASTA + metadata |
| `snakemake-minimal` | Only for the Snakemake route | Runs `workflow/Snakefile` over many samples (install from bioconda; the pip package is `snakemake`) |
| `nextflow` | Only for the Nextflow route | Runs `workflow/main.nf` over many samples (bioconda; needs Java 17+, provided by `openjdk`) |
| `r-base` (env `ngs-r`) | For step 9b | Runs `R/variant_stats.R` (statistics shown in the dashboard) |

If not using the conda command above:

```bash
# polars    -> VCF parsing/filtering into the final TSV report
# tqdm      -> download progress bars
# requests  -> ENA API + FASTQ download over HTTPS
# biopython -> NCBI Entrez access (reference genome + metadata)
pip install polars tqdm requests biopython
```

## 3. SnpEff database

SnpEff can only annotate a variant if its database describes the **same reference sequence** (same chromosome name and coordinates). The script therefore picks the database from the `SNPEFF_DB` table, keyed by reference accession:

| Reference | SnpEff database | Status |
|-----------|-----------------|--------|
| `NC_000913.3` (E. coli K-12) | `Escherichia_coli_K_12_substr_MG1655` | Mapped by default |
| Any other (e.g. `NC_002549.1` Ebola) | — | Annotation is **skipped** with a warning; the report's Effect/Impact/Gene columns stay empty |

Download the E. coli database once before running. If it is missing, annotation fails, the log tells you which download command to run, and the report is built from the raw VCF:

```bash
# Download the prebuilt E. coli K-12 MG1655 annotation database (genes, transcripts, codons)
# -v = verbose output
snpEff download -v Escherichia_coli_K_12_substr_MG1655

# Download the prebuilt Ebola Zaire virus annotation database (genes, transcripts, codons)
# Database: ebola_zaire (based on KJ660346.1 / Makona strain)
# -v = verbose output
snpEff download -v ebola_zaire
```

For other organisms, list available databases and add an `"<reference accession>": "<database>"` entry to `SNPEFF_DB` in the script:

```bash
# Find a candidate database for your organism
snpEff databases | grep -i "<organism>"

# Check that the database's chromosome name matches your reference FASTA header
snpEff download -v <database>
snpEff dump <database> | head
grep ">" ngs-variant-calling/<run>/01_reference/<accession>.fasta
```

Then add the pair to `SNPEFF_DB` in `run_interactive_pipeline.py`:

```python
SNPEFF_DB = {
    "NC_000913.3": "Escherichia_coli_K_12_substr_MG1655",
    "<accession>": "<database>",
}
```

## 4. Verify installation

```bash
# Check that every external tool the pipeline calls is on PATH
for t in fastp bwa samtools bcftools bgzip tabix snpEff java; do
    command -v "$t" >/dev/null && echo "OK       $t" || echo "MISSING  $t"
done
# Check that the Python libraries import (Bio = biopython)
python3 -c "import polars, tqdm, requests, Bio; print('Python packages OK')"
# Snakemake route only: version check, then a dry run (lists the jobs, runs nothing)
snakemake --version
cd ngs-variant-calling && snakemake -s workflow/Snakefile -n && cd ..
# Nextflow route only
nextflow -version
# R for step 9b (prints the Rscript the pipeline will use)
python3 -c "import sys; sys.path.insert(0,'ngs-variant-calling'); import run_interactive_pipeline as m; print(m.find_rscript() or 'MISSING Rscript')"
```

## 5. System requirements

- A web browser to open the results dashboard (`07_dashboard/dashboard.html`). It needs no internet and no extra packages: `dashboard.py` uses only the Python standard library plus `samtools`/`bcftools`.
- Optional: Power BI Desktop or Excel. The dashboard already shows its CSV tables in the browser; these are only needed for your own analysis.

- Linux (Ubuntu or GitHub Codespaces recommended)
- Internet access (ENA / NCBI downloads)
- ~2–4 GB free disk space per run

## 6. How this file was created

This file was generated from the terminal with a bash **heredoc**, which writes multi-line text straight into a file without opening an editor:

````bash
cd /workspaces/Bioinformatic_workflow

cat > Prerequisites.md <<'EOF'
# Prerequisites

Everything needed to run `ngs-variant-calling/run_interactive_pipeline.py`.

## 1. Command-line tools (via Conda/Mamba, bioconda channel)
...
(rest of the Markdown content)
...
EOF

# Confirm the file was written
ls -l Prerequisites.md
cat Prerequisites.md
````

How it works:

- `cat > Prerequisites.md` sends everything that follows into `Prerequisites.md`, creating the file or **overwriting** it.
- `<<'EOF'` starts the heredoc, and the text continues until a line that contains only `EOF`.
- Putting `'EOF'` in quotes stops the shell from expanding `$VARIABLES`, `` `backticks` `` or `$(commands)` inside the text, so lines such as `source "$HOME/miniforge3/bin/activate"` are written literally.
- To **add** text to an existing file instead of overwriting it, use `>>`. This section was appended that way:

  ```bash
  cat >> Prerequisites.md <<'OUTER'
  ## 6. How this file was created
  ...
  OUTER
  ```

  A different end marker (`OUTER`) was used here because the appended text itself contains a line with only `EOF`.

Later updates (comments, troubleshooting, the `SNPEFF_DB` notes) were made by editing the file in place rather than by re-running the heredoc, because `cat >` would have overwritten the whole file.
