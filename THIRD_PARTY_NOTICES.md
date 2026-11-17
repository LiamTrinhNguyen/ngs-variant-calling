# Third-party software and data

The code and documentation in this repository are licensed under the MIT License (see [`LICENSE`](LICENSE)). That covers:
- the pipeline code: `run_interactive_pipeline.py`, `dashboard.py`, `R/variant_stats.R`, the Snakemake and Nextflow workflows and `run_nextflow.sh`;
- the guides, diagrams and pages: the `.md` files, `.svg` files, `index.html` and `showcase.html`.

The pipeline also relies on other people's software and public data. None of it is copied into this repository, and the MIT License does not cover it. Each part keeps its own license, as listed below.

This file is a summary to help you use the pipeline responsibly. It is not legal advice.

## Software the pipeline runs

These are installed separately into conda environments (see `Prerequisites.md`), not shipped here. The command-line tools run as separate programs. The Python packages are imported by the scripts.

| Software | Version tested | License | Used for |
|---|---|---|---|
| [BWA](https://github.com/lh3/bwa) | 0.7.19 | GPL-3.0-or-later | read alignment (`bwa index`, `bwa mem`) |
| [SAMtools](https://www.htslib.org/) | 1.25 | MIT | sorting, indexing, `faidx`, `dict` |
| [HTSlib](https://www.htslib.org/) (`bgzip`, `tabix`) | 1.25 | MIT | compressing and indexing VCFs |
| [BCFtools](https://www.htslib.org/) | 1.25 | GPL-3.0-or-later (the conda build includes GSL) | variant calling, normalising, flagging |
| [fastp](https://github.com/OpenGene/fastp) | 1.4.0 | MIT | read QC and trimming |
| [SnpEff](https://pcingola.github.io/SnpEff/) | 5.4.0c | LGPL-3.0 | variant annotation (optional) |
| [OpenJDK](https://openjdk.org/) | 25.0.2 | GPL-2.0-or-later with Classpath exception | runs SnpEff and Nextflow |
| [Snakemake](https://snakemake.github.io/) (`snakemake-minimal`) | 9.27.1 | MIT | multi-sample workflow |
| [Nextflow](https://www.nextflow.io/) | 26.10.6 | Apache-2.0 | multi-sample workflow |
| [R](https://www.r-project.org/) (`r-base`, env `ngs-r`) | 4.6.2 | GPL-2.0-or-later | statistics (step 9b), base R only |
| [Kraken2](https://github.com/DerrickWood/kraken2) (env `ngs-kraken`, optional) | 2.17.2 | GPL-3.0-or-later | read classification (step 5d) |
| [Bracken](https://github.com/jenniferlu717/Bracken) (env `ngs-kraken`, optional) | 3.1 | GPL-3.0-or-later | species abundance re-estimation (step 5d) |
| [Python](https://www.python.org/) | 3.11 | PSF License | runs the scripts |
| [Polars](https://pola.rs/) | 2.1.0 | MIT | tables and reports |
| [Biopython](https://biopython.org/) | 1.88 | Biopython License Agreement / BSD-3-Clause | downloading references from NCBI |
| [Requests](https://requests.readthedocs.io/) | 2.34.2 | Apache-2.0 | downloading reads from ENA |
| [tqdm](https://tqdm.github.io/) | 4.70.1 | MPL-2.0 and MIT | progress bars |
| [PyYAML](https://pyyaml.org/) | 6.0.3 | MIT | reading `config.yaml` |
| [Miniforge](https://github.com/conda-forge/miniforge) | – | BSD-3-Clause | installs conda/mamba |

The licenses come from the conda package metadata of the versions listed (`~/miniforge3/envs/<env>/conda-meta/*.json`). Newer versions may change them. Check that metadata, or each project's own license file, before redistributing.

**About the GPL tools (BWA, BCFtools, R, Kraken2, Bracken):** this repository only calls them as separate programs. It does not contain, modify or link to their code, so their licenses don't change how you can use the MIT-licensed code here. If you **redistribute** the tools themselves, for example in a Docker image or a bundled installer, you must follow their licenses, which includes making their source code available.

## Data

- **Sequencing reads** are downloaded at run time from the [European Nucleotide Archive (ENA)](https://www.ebi.ac.uk/ena/browser/), which mirrors NCBI SRA and DDBJ through the International Nucleotide Sequence Database Collaboration (INSDC).
- **Reference genomes** are downloaded at run time from [NCBI Nucleotide / RefSeq](https://www.ncbi.nlm.nih.gov/nuccore).
- **SnpEff databases** are downloaded by `snpEff download` from the SnpEff project.
- **Kraken2 databases** (optional, step 5d) are downloaded by you, e.g. the prebuilt indexes at https://benlangmead.github.io/aws-indexes/k2. They are built from NCBI RefSeq and other public sources; check the page for each database's sources and terms.
- **A host genome** for host removal (optional, step 5b) is supplied by you. `Prerequisites.md` suggests the human T2T-CHM13 assembly from the Human Pangenome Reference Consortium's public bucket; check the consortium's data-use terms before redistributing it.
- **The reference screen** (step 5c, and option `0` at question 4 of the interactive script) downloads the built-in references and two human decoy sequences (the rDNA repeat `U13369.1` and the mitochondrial genome `NC_012920.1`) from NCBI once into `~/.cache/ngs_pipeline/reference_screen/`, and for option `0` on an SRA run streams the first 50,000 reads from ENA. Neither is stored in this repository.

INSDC and NCBI make public sequence data available without restrictions on use. However, NCBI notes that original submitters may claim patent, copyright or other rights in the data they submitted. You are responsible for checking the terms for any dataset you use.

The example run outputs in this repository (`ngs-variant-calling/NGS_run/` and the other run folders) are derived from the public run **SRR1553607** (Zaire ebolavirus) and the reference **NC_002549.1**. If you publish results based on them, cite the original study listed on the run's ENA or SRA page.

## Screenshots

The images in `image/` are screenshots of Claude Code, taken for the setup guide `integrate_claude.md`. Claude Code and its interface belong to Anthropic and are not covered by this repository's MIT License.

## Citing the tools

If you publish results from this pipeline, please cite the tools that produced them, at least BWA, SAMtools/BCFtools, fastp and, if used, SnpEff, Snakemake or Nextflow, and R. Each project's website lists its preferred citation.
