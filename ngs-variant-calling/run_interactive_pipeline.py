#!/usr/bin/env python3
"""
NGS_Pipeline_v2.6 – Full automated variant calling pipeline
Author : Liam TrinhNguyen
Email  : LiamTrinhNguyen@gmail.com

Three ways to run it:
  1. Interactive, one sample:   python3 run_interactive_pipeline.py
  2. Snakemake, many samples:   snakemake -s workflow/Snakefile -c 4
  3. Nextflow, many samples:    nextflow run workflow/main.nf
     (2 and 3 both read their samples from config/config.yaml)

The Snakefile and workflow/main.nf import this module and call the step
functions below (get_sra_metadata_and_fastq, download_reference,
index_reference, run_fastp, align_and_sort, call_variants, annotate_variants,
filter_and_report, run_r_stats) and the SNPEFF_DB table. If you change a
function's signature or the names of the files it writes, update the matching
Snakemake rule AND Nextflow process as well.

Step 9b (statistics) is the R script R/variant_stats.R, run by run_r_stats().

Step 10 (the results dashboard) lives in dashboard.py; both routes call
dashboard.build_dashboard() at the end.
"""

__author__  = "Liam TrinhNguyen"
__email__   = "LiamTrinhNguyen@gmail.com"
__version__ = "NGS_Pipeline_v2.6"

import os
import json  # saves the calling settings next to each VCF
import re   # used to match R1/R2 file names
import sys
import logging
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Tuple

import polars as pl

# ---------------------------------------------------------------------------
# Optional libraries
#   WHY: the script still starts if one is missing, and only the step that
#        needs it fails, with a message saying what to install.
# ---------------------------------------------------------------------------
# tqdm – progress bars for downloads (installed automatically if missing)
try:
    from tqdm import tqdm
except ImportError:
    print("Installing tqdm...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "tqdm", "--quiet"])
    from tqdm import tqdm

# requests – downloads reads and metadata from ENA (SRA accessions)
try:
    import requests
except ImportError:
    requests = None

# Biopython – downloads the reference genome from NCBI
try:
    from Bio import Entrez, SeqIO
    Entrez.email = "LiamTrinhNguyen@gmail.com"
except ImportError:
    Entrez = None
    SeqIO = None

# ---------------------------------------------------------------------------
# Logging
#   WHAT: prints time | level | message for every step.
#   WHY:  uses its own logger instead of logging.basicConfig(), because
#         basicConfig changes Python's shared logger, and Snakemake then
#         printed every line twice.
# ---------------------------------------------------------------------------
log = logging.getLogger("NGS_Pipeline")
if not log.handlers:
    _handler = logging.StreamHandler()
    _handler.setFormatter(logging.Formatter(
        "%(asctime)s | %(levelname)-8s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    ))
    log.addHandler(_handler)
log.setLevel(logging.INFO)
log.propagate = False

# ---------------------------------------------------------------------------
# STEP 1 · SETUP – create the output folders
#   WHAT: 00_raw … 07_dashboard + logs under the output folder.
#   WHY:  each step writes to its own numbered folder, so results are easy
#         to find and separate runs don't mix.
# ---------------------------------------------------------------------------
def setup_directories(base_dir: Path) -> dict:
    """Create the numbered output folders and return their paths."""
    dirs = {
        "base":      base_dir,
        "raw":       base_dir / "00_raw",
        "ref":       base_dir / "01_reference",
        "qc":        base_dir / "02_qc",
        "aligned":   base_dir / "03_aligned",
        "variants":  base_dir / "04_variants",
        "annotated": base_dir / "05_annotated",
        "reports":   base_dir / "06_reports",
        "dashboard": base_dir / "07_dashboard",
        "logs":      base_dir / "logs",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    log.info(f"Project directories created under {base_dir}")
    return dirs

# ---------------------------------------------------------------------------
# Helper · run an external tool (bwa, samtools, bcftools, …)
#   WHY: one place that logs the command, saves the tool's output to a log
#        file, and stops with a clear error if the tool fails.
# ---------------------------------------------------------------------------
def run_cmd(cmd: List[str], log_file: Optional[Path] = None, check: bool = True) -> subprocess.CompletedProcess:
    """Execute a shell command and optionally capture stdout/stderr to a log."""
    log.info("CMD: " + " ".join(str(c) for c in cmd))
    if log_file:
        with open(log_file, "w") as fh:
            result = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, text=True)
    else:
        result = subprocess.run(cmd, capture_output=True, text=True)

    if check and result.returncode != 0:
        err = result.stderr or result.stdout or "No error output"
        raise RuntimeError(f"Command failed (exit {result.returncode}):\n{err}")
    return result

# ---------------------------------------------------------------------------
# STEP 4 · DATA INPUT – get the reads
#   WHAT: an SRA run (SRR/ERR/DRR) is downloaded from ENA, and its metadata
#         (organism, platform, read count) is printed.
#   WHY:  the organism name is used later to pick the reference genome
#         automatically, and the metadata lets you check it's the right sample.
# ---------------------------------------------------------------------------
def process_accession(accession_id: str, dirs: dict) -> Tuple[str, List[Path]]:
    """
    Detect SRA vs Nuccore accession and fetch data.
    Returns (organism_name, list_of_fastq_paths)
    """
    accession_id = accession_id.strip()
    if accession_id.startswith(("SRR", "ERR", "DRR", "SAM", "PRN", "SRS", "ERS", "DRS")):
        log.info(f"[{accession_id}] Detected SRA / Read Archive Accession")
        return get_sra_metadata_and_fastq(accession_id, dirs["raw"])
    else:
        log.info(f"[{accession_id}] Detected Nuccore / Reference Accession")
        get_ncbi_metadata(accession_id)
        fasta = download_ncbi_fasta(accession_id, dirs["ref"] / f"{accession_id}.fasta")
        return "unknown", [fasta]

def get_ncbi_metadata(accession_id: str) -> None:
    if Entrez is None:
        log.warning("Biopython not available – skipping detailed metadata")
        return
    log.info(f"Fetching metadata for {accession_id}")
    try:
        handle = Entrez.esummary(db="nuccore", id=accession_id, retmode="json")
        record = Entrez.read(handle)
        handle.close()
        uid = record["uids"][0]
        doc = record[uid]

        organism = doc.get("Organism", "N/A")
        scientific_name = doc.get("ScientificName", "N/A")
        tax_id = doc.get("TaxId", "N/A")
        description = doc.get("Title", "N/A")
        seq_length = doc.get("Length", "N/A")

        seq_handle = Entrez.efetch(db="nuccore", id=accession_id, rettype="fasta", retmode="text")
        fasta_record = SeqIO.read(seq_handle, "fasta")
        seq_handle.close()
        sequence = str(fasta_record.seq).upper()
        total_len = len(sequence)
        a = sequence.count("A")
        c = sequence.count("C")
        g = sequence.count("G")
        t = sequence.count("T")
        other = total_len - (a + c + g + t)
        gc = ((g + c) / total_len * 100) if total_len else 0.0

        print("\n" + "=" * 55)
        print(f"METRICS & METADATA FOR: {accession_id}")
        print("=" * 55)
        print(f"Description     : {description}")
        print(f"Organism        : {organism}")
        print(f"Scientific Name : {scientific_name}")
        print(f"Taxonomy ID     : {tax_id}")
        print(f"Sequence Length : {seq_length:,} bp")
        print(f"GC Content      : {gc:.2f}%")
        print(f"Base Counts     : A:{a:,}  C:{c:,}  G:{g:,}  T:{t:,}  Other:{other:,}")
        print("=" * 55 + "\n")
    except Exception as e:
        log.error(f"Metadata fetch failed: {e}")

def download_ncbi_fasta(accession_id: str, output_path: Path) -> Path:
    if Entrez is None:
        raise RuntimeError(
            "Biopython is required to download NCBI sequences.\n"
            "Install it with:\n"
            "    pip install biopython"
        )
    log.info(f"Downloading FASTA → {output_path}")
    handle = Entrez.efetch(db="nuccore", id=accession_id, rettype="fasta", retmode="text")
    data = handle.read()
    handle.close()
    output_path.write_text(data)
    log.info(f"Saved {output_path}")
    return output_path

def get_sra_metadata_and_fastq(run_accession: str, out_dir: Path) -> Tuple[str, List[Path]]:
    """Fetch ENA metadata + download FASTQ via HTTPS. Returns (organism, [fastq paths])."""
    if requests is None:
        raise RuntimeError(
            "The 'requests' library is required for ENA downloads.\n"
            "Install it with:\n"
            "    pip install requests"
        )

    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?"
        f"accession={run_accession}&result=read_run&"
        "fields=study_accession,sample_accession,experiment_accession,"
        "scientific_name,instrument_platform,instrument_model,"
        "library_layout,library_strategy,read_count,base_count,fastq_ftp"
    )
    log.info(f"Querying ENA for {run_accession}")
    resp = requests.get(url, timeout=90)
    if resp.status_code != 200:
        raise RuntimeError(f"ENA API failed: HTTP {resp.status_code}")

    lines = resp.text.strip().splitlines()
    if len(lines) < 2:
        raise RuntimeError(f"No metadata returned for {run_accession}")

    header = lines[0].split("\t")
    values = lines[1].split("\t")
    data = dict(zip(header, values))

    organism = data.get("scientific_name", "unknown")
    platform = data.get("instrument_platform", "N/A")
    model    = data.get("instrument_model", "N/A")
    strategy = data.get("library_strategy", "N/A")
    layout   = data.get("library_layout", "N/A")
    reads    = data.get("read_count", "N/A")
    bases    = data.get("base_count", "N/A")

    print("\n" + "=" * 60)
    print(f"SRA METADATA & STATISTICS FOR: {run_accession}")
    print("=" * 60)
    print(f"Scientific Name    : {organism}")
    print(f"Sequencing Platform: {platform} ({model})")
    print(f"Library Strategy   : {strategy} ({layout})")
    if str(reads).isdigit():
        print(f"Total Reads        : {int(reads):,}")
    if str(bases).isdigit():
        print(f"Total Bases        : {int(bases):,}")
    print("=" * 60 + "\n")

    fastq_paths = []
    for link in data.get("fastq_ftp", "").split(";"):
        if not link.strip():
            continue
        clean = link.replace("ftp://", "").strip()
        file_url = f"https://{clean}"
        fname = out_dir / Path(file_url).name
        log.info(f"Downloading {fname.name} …")

        with requests.get(file_url, stream=True, timeout=600) as r:
            r.raise_for_status()
            total = int(r.headers.get("content-length", 0))
            with open(fname, "wb") as f, tqdm(
                total=total, unit="B", unit_scale=True, desc=fname.name, ncols=80
            ) as bar:
                for chunk in r.iter_content(chunk_size=1 << 20):
                    if chunk:
                        f.write(chunk)
                        bar.update(len(chunk))
        fastq_paths.append(fname)
        log.info(f"Saved {fname}")

    if not fastq_paths:
        raise RuntimeError(f"No FASTQ files found for {run_accession}")

    # Pick R1/R2 by file name (_1 / _2), not by list position.
    # WHY: ENA can list 3 files (SRR.fastq.gz, SRR_1, SRR_2). Taking the first
    #      two would pair the leftover unpaired file with R1 and give wrong
    #      alignments. No _1/_2 pair → treat the run as single-end.
    mate1 = next((p for p in fastq_paths if p.name.endswith("_1.fastq.gz")), None)
    mate2 = next((p for p in fastq_paths if p.name.endswith("_2.fastq.gz")), None)
    if mate1 and mate2:
        return organism, [mate1, mate2]
    return organism, [fastq_paths[0]]

# ---------------------------------------------------------------------------
# STEP 2 · REFERENCE – download the reference genome
#   WHAT: the accession you gave, or one looked up from the organism name in
#         REF_MAP, downloaded from NCBI as FASTA into 01_reference/.
#   WHY:  variants are the differences between your reads and a known
#         genome, so a reference is needed to compare against.
# ---------------------------------------------------------------------------
def download_reference(organism: str, ref_dir: Path, force_accession: str = None) -> Path:
    """
    Download a suitable reference genome.
    You can force a specific accession with force_accession.
    """
    if force_accession:
        acc = force_accession.strip()
        log.info(f"Using forced reference accession: {acc}")
    else:
        # Organism name → RefSeq accession. Matching ignores case and works on
        # part of the name ("Zaire ebolavirus" matches "Zaire ebolavirus Makona").
        # To support a new organism, add a line here.
        REF_MAP = {
            # Bacteria
            "Escherichia coli": "NC_000913.3",
            "Escherichia coli K-12": "NC_000913.3",
            "Salmonella enterica": "NC_003198.1",
            "Staphylococcus aureus": "NC_007795.1",
            "Mycobacterium tuberculosis": "NC_000962.3",

            # Viruses – Ebola family
            "Zaire ebolavirus": "NC_002549.1",          # Mayinga (most used)
            "Ebola virus": "NC_002549.1",
            "Zaire ebolavirus Mayinga": "NC_002549.1",
            "Sudan ebolavirus": "NC_006432.1",
            "Bundibugyo ebolavirus": "NC_014373.1",
            "Tai Forest ebolavirus": "NC_014372.1",
            "Reston ebolavirus": "NC_004161.1",
            "Marburg marburgvirus": "NC_001608.3",

            # Other viruses
            "SARS-CoV-2": "NC_045512.2",
            "Influenza A virus": "NC_002016.1",

            # Eukaryotes (examples)
            "Homo sapiens": "GCF_000001405.40",
            "Saccharomyces cerevisiae": "GCF_000146045.2",
        }

        key = next((k for k in REF_MAP if k.lower() in organism.lower()), None)

        if key is None:
            raise RuntimeError(
                f"\nNo pre-mapped reference genome found for organism: '{organism}'.\n"
                f"Please either:\n"
                f"  1. Add the organism → accession mapping to REF_MAP, or\n"
                f"  2. Force a reference with force_accession='NC_xxxxxx.x'\n"
            )
        acc = REF_MAP[key]
        log.info(f"Selected reference for '{organism}': {acc}  ({key})")

    fasta = ref_dir / f"{acc}.fasta"
    # WHY: skip the download if the file is already there (saves time on re-runs)
    if not fasta.exists():
        download_ncbi_fasta(acc, fasta)
    else:
        log.info(f"Reference already present: {fasta}")

    return fasta

# ---------------------------------------------------------------------------
# STEP 3 · INDEXING – prepare the reference for fast lookup
#   WHAT: bwa index (.bwt .pac .sa …), samtools faidx (.fai), samtools dict (.dict).
#   WHY:  like the index of a book: aligners and variant callers can jump to
#         any position instead of reading the whole genome every time.
# ---------------------------------------------------------------------------
def index_reference(fasta: Path, log_dir: Path) -> None:
    # bwa index – needed by BWA-MEM to place reads on the genome
    log.info(f"Building BWA index for {fasta.name}")
    run_cmd(["bwa", "index", str(fasta)], log_dir / "bwa_index.log")

    # samtools faidx – lets bcftools read the reference base at any position
    log.info("Building samtools faidx")
    run_cmd(["samtools", "faidx", str(fasta)], log_dir / "samtools_faidx.log")

    # samtools dict – lists chromosome names and lengths (used by many tools)
    dict_file = fasta.with_suffix(".dict")
    if not dict_file.exists():
        run_cmd(
            ["samtools", "dict", str(fasta), "-o", str(dict_file)],
            log_dir / "samtools_dict.log",
        )

# ---------------------------------------------------------------------------
# STEP 5 · QC & TRIMMING – fastp
#   WHAT: removes adapters and low-quality bases, drops very short reads,
#         and writes an HTML/JSON quality report to 02_qc/.
#   WHY:  adapters and sequencing errors would otherwise look like
#         mutations, so they're cleaned up before alignment.
# ---------------------------------------------------------------------------
def run_fastp(
    r1: Path,
    r2: Optional[Path],
    out_dir: Path,
    sample: str,
    threads: int = 4,
) -> Tuple[Path, Optional[Path]]:
    out_r1 = out_dir / f"{sample}_R1.trimmed.fastq.gz"
    out_r2 = out_dir / f"{sample}_R2.trimmed.fastq.gz" if r2 else None
    html   = out_dir / f"{sample}.fastp.html"
    json   = out_dir / f"{sample}.fastp.json"

    cmd = [
        "fastp",
        "-i", str(r1),
        "-o", str(out_r1),
        "-h", str(html),
        "-j", str(json),
        "--correction",                     # fix mismatches where R1 and R2 overlap
        "--thread", str(threads),
        "--qualified_quality_phred", "20",  # a base is "good" at Q20 (99% accurate)
        "--length_required", "30",          # drop reads < 30 bp: too short to place uniquely
    ]
    # Find adapters automatically, so you don't need to know the kit used
    if r2:
        cmd.extend(["-I", str(r2), "-O", str(out_r2), "--detect_adapter_for_pe"])
    else:
        cmd.append("--detect_adapter_for_se")

    run_cmd(cmd, out_dir / f"{sample}.fastp.log")
    return out_r1, out_r2

# ---------------------------------------------------------------------------
# STEP 6 · ALIGNMENT – BWA-MEM + samtools
#   WHAT: finds where each read belongs on the reference, then sorts and
#         indexes the result as 03_aligned/<sample>.sorted.bam.
#   WHY:  reads must be lined up on the genome before you can see which
#         positions differ from the reference.
# ---------------------------------------------------------------------------
def align_and_sort(
    ref: Path,
    r1: Path,
    r2: Optional[Path],
    out_dir: Path,
    sample: str,
    threads: int = 8,
) -> Path:
    sam = out_dir / f"{sample}.sam"
    bam = out_dir / f"{sample}.sorted.bam"

    bwa_cmd = [
        "bwa", "mem",
        "-t", str(threads),
        # Read group: tags every read with the sample name.
        # WHY: the VCF uses SM as the sample column name.
        "-R", f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA",
        str(ref), str(r1)
    ]
    if r2:
        bwa_cmd.append(str(r2))

    log.info("Running BWA-MEM")
    with open(sam, "w") as fh:
        subprocess.run(bwa_cmd, stdout=fh, check=True)

    # Sort by position + index (.bai)
    # WHY: variant callers walk the genome in order and need to jump to regions
    log.info("Sorting & indexing BAM")
    run_cmd(["samtools", "sort", "-@", str(threads), "-o", str(bam), str(sam)])
    run_cmd(["samtools", "index", str(bam)])
    # WHY: the SAM is a large uncompressed copy of the BAM; delete it to save disk
    sam.unlink(missing_ok=True)
    return bam

# ---------------------------------------------------------------------------
# STEP 7 · VARIANT CALLING – bcftools
#   WHAT: mpileup counts the bases seen at each position; call decides where
#         they differ from the reference → 04_variants/<sample>.raw.vcf.gz.
#   WHY:  this is the step that actually finds the mutations.
# ---------------------------------------------------------------------------
# Settings for calling (step 7) and filtering (step 9), kept in one place.
# Snakemake can override the three filters from config.yaml.
MIN_MAPQ  = 20      # ignore reads with mapping quality < 20
MIN_BASEQ = 20      # ignore bases with quality < 20 (99% accurate)
MAX_DEPTH = 10000   # reads looked at per position
MIN_QUAL  = 30.0    # keep: caller confidence  QUAL >= 30
MIN_DP    = 15      # keep: good-quality reads DP   >= 15
MIN_AF    = 0.8     # keep: share of reads with the variant  AF >= 0.8


def call_variants(ref: Path, bam: Path, out_dir: Path, sample: str) -> Path:
    vcf = out_dir / f"{sample}.raw.vcf.gz"
    log.info("Calling variants with bcftools mpileup | call | norm")

    # mpileup – count the bases at each position
    #   -q 20 / -Q 20 : only count reads with MAPQ >= 20 and bases with quality >= 20.
    #     WHY: a read that fits equally well in several places (MAPQ 0), or a
    #          blurry base, is weak evidence and can create false variants.
    #   -d 10000 : look at up to 10,000 reads per position.
    #     WHY: the default (250) silently ignores the rest; viral and amplicon
    #          data often has thousands of reads per position.
    #   -a FORMAT/AD,FORMAT/DP : record reads per allele (AD) and good-quality depth (DP).
    #     WHY: needed to compute the allele frequency (AF) used in the filter.
    # call – decide where the sample differs from the reference
    #   --ploidy 1 : haploid. WHY: bacteria and viruses have one genome copy, so a
    #                site is either reference or variant; diploid calling gave
    #                fake "0/1" half-calls. Remove this flag for human, yeast, etc.
    #   -m : multiallelic caller (recommended model)   -v : output variant sites only
    # norm – write every variant in one standard form
    #   -f ref : left-align indels.  -m -any : split multi-allelic sites into one line each.
    #     WHY: the same indel can be written several ways; without this, the same
    #          variant can look different between samples or databases.
    cmd = (
        f"bcftools mpileup -f {ref} -q {MIN_MAPQ} -Q {MIN_BASEQ} -d {MAX_DEPTH} "
        f"-a FORMAT/AD,FORMAT/DP -Ou {bam} | "
        f"bcftools call --ploidy 1 -mv -Ou | "
        f"bcftools norm -f {ref} -m -any -Oz -o {vcf}"
    )
    run_cmd(["bash", "-c", cmd], out_dir / f"{sample}.bcftools.log")
    # WHY: the index (.csi) lets later tools jump to any region of the VCF
    run_cmd(["bcftools", "index", "-f", str(vcf)])

    # Save the settings used, next to the VCF.
    # WHY: bcftools mpileup does not write its options into the VCF header, so
    #      without this file nobody (including the dashboard) could tell later
    #      which quality cut-offs this VCF was made with.
    (out_dir / f"{sample}.calling_settings.json").write_text(json.dumps({
        "min_mapq": MIN_MAPQ, "min_baseq": MIN_BASEQ, "max_depth": MAX_DEPTH,
        "haploid": True, "normalized": True, "allele_depths": True,
        "command": cmd,
    }, indent=2))
    return vcf

# ---------------------------------------------------------------------------
# STEP 8 · ANNOTATION – SnpEff
#   WHAT: adds what each variant does: gene, effect (e.g. missense),
#         impact (HIGH/MODERATE/LOW/MODIFIER) → 05_annotated/.
#   WHY:  a position alone says little; the annotation tells you whether a
#         gene or protein is changed.
# ---------------------------------------------------------------------------
# Reference accession → SnpEff database name.
# WHY: SnpEff finds variants by chromosome name, so the database must
#      describe the SAME reference. A wrong database annotates nothing.
#      Check before adding one:  snpEff databases | grep -i <organism>
SNPEFF_DB = {
    "NC_000913.3": "Escherichia_coli_K_12_substr_MG1655",   # E. coli K-12 MG1655
}

def annotate_variants(vcf: Path, genome: str, out_dir: Path, sample: str) -> Path:
    """Annotate a VCF with the SnpEff database `genome` (see SNPEFF_DB)."""
    ann_vcf = out_dir / f"{sample}.ann.vcf.gz"
    log.info(f"Annotating with SnpEff ({genome})")

    # SnpEff writes plain text, so write it uncompressed first
    uncompressed = ann_vcf.with_suffix("")
    cmd = ["snpEff", "-v", genome, str(vcf)]
    try:
        with open(uncompressed, "w") as fh:
            subprocess.run(cmd, stdout=fh, check=True)
    except Exception:
        # WHY: a failed run leaves an empty file that looks like a result;
        #      delete it so the failure is obvious
        uncompressed.unlink(missing_ok=True)
        raise

    # WHY: compress (bgzip) and index (tabix) to save space and allow quick lookups
    run_cmd(["bgzip", "-f", str(uncompressed)])
    run_cmd(["tabix", "-p", "vcf", str(ann_vcf)])
    return ann_vcf

# ---------------------------------------------------------------------------
# STEP 9 · FILTER & REPORT – bcftools filter + Polars
#   WHAT: FLAGS every call with the reason it fails (it doesn't delete it),
#         then writes two tables:
#           06_reports/<sample>.flagged.vcf.gz             all calls, FILTER column set
#           06_reports/<sample>.all_calls.tsv              all calls + FILTER reason
#           06_reports/<sample>.high_confidence_variants.tsv   only FILTER = PASS
#   WHY:  flagging keeps an audit trail. You can see what was excluded and why,
#         which labs need to justify a result. The PASS table holds only the
#         variants worth looking at.
# ---------------------------------------------------------------------------
REPORT_COLS = ["CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AF",
               "Effect", "Impact", "Gene", "HGVS_c", "HGVS_p"]


def flag_calls(vcf: Path, out_vcf: Path, min_qual: float, min_dp: int, min_af: float) -> Path:
    """Set the VCF FILTER column: PASS, or the reasons (LowQual;LowDepth;LowAF)."""
    # Three soft filters in a row: -s names the reason, -m + adds it to earlier ones.
    #   LowQual  : QUAL < min_qual          – caller not confident enough
    #   LowDepth : FORMAT/DP < min_dp       – too few GOOD-quality reads
    #              WHY FORMAT/DP, not INFO/DP: INFO/DP counts every read, even
    #              the low-quality ones that -q/-Q told mpileup to ignore.
    #   LowAF    : AD[alt] / DP < min_af    – too few reads show the variant
    #              WHY: in a haploid genome a real variant is in (nearly) all
    #              reads. 50–80% suggests a mixed sample, contamination or a
    #              mapping error, not a clean call.
    cmd = (
        f"bcftools filter -s LowQual -e 'QUAL<{min_qual}' {vcf} -Ou | "
        f"bcftools filter -m + -s LowDepth -e 'FMT/DP<{min_dp}' -Ou | "
        f"bcftools filter -m + -s LowAF -e 'FMT/AD[0:1]/FMT/DP<{min_af}' -Oz -o {out_vcf}"
    )
    run_cmd(["bash", "-c", cmd])
    run_cmd(["bcftools", "index", "-f", str(out_vcf)])
    return out_vcf


def filter_and_report(
    vcf: Path,
    reports_dir: Path,
    sample: str,
    min_qual: float = MIN_QUAL,
    min_dp: int = MIN_DP,
    min_af: float = MIN_AF,
) -> Path:
    """
    Flag every call, then write all_calls.tsv (with FILTER) and
    high_confidence_variants.tsv (PASS only). Works on annotated or raw VCFs.
    Returns the path of the high-confidence report.
    """
    log.info(f"Flagging calls: QUAL >= {min_qual}, DP >= {min_dp}, AF >= {min_af}")
    flagged = flag_calls(vcf, reports_dir / f"{sample}.flagged.vcf.gz", min_qual, min_dp, min_af)

    tmp_tsv = reports_dir / f"{sample}.tmp.tsv"
    # Pull the columns out of the VCF. [%DP] / [%AD] are the sample's good-quality
    # depth and reads per allele (ref,alt). Try with SnpEff's ANN field first; a raw
    # VCF has no ANN, so bcftools fails and we retry without it.
    base_fmt = "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%FILTER\t[%DP]\t[%AD]"
    try:
        run_cmd(["bcftools", "query", "-f", base_fmt + "\t%INFO/ANN\n", str(flagged), "-o", str(tmp_tsv)])
        has_ann = True
    except RuntimeError:
        log.warning("No INFO/ANN field found – treating as raw (unannotated) VCF")
        run_cmd(["bcftools", "query", "-f", base_fmt + "\n", str(flagged), "-o", str(tmp_tsv)])
        has_ann = False

    cols = ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "DP", "AD"] + (["ANN"] if has_ann else [])
    df = pl.read_csv(
        tmp_tsv, separator="\t", has_header=False, new_columns=cols,
        schema_overrides={c: pl.Utf8 for c in ["CHROM", "REF", "ALT", "FILTER", "AD"]},
        ignore_errors=True,
    )

    # AF = reads with the variant / good-quality reads at the site.
    # AD is "ref,alt" (one ALT per line after bcftools norm -m -any).
    alt_reads = pl.col("AD").str.split(",").list.get(1, null_on_oob=True).cast(pl.Int64, strict=False)
    dp = pl.col("DP").cast(pl.Int64, strict=False)
    df = df.with_columns(
        pl.when(dp > 0).then((alt_reads / dp).round(3)).otherwise(None).alias("AF")
    )

    if has_ann:
        # Split the SnpEff ANN text into columns. Its fields are separated by "|":
        #   0 Allele | 1 Effect | 2 Impact | 3 Gene | … | 9 HGVS.c | 10 HGVS.p
        # A variant can have several annotations (separated by ",").
        # WHY the first one: SnpEff lists the most severe annotation first.
        first_ann = pl.col("ANN").str.split(",").list.get(0, null_on_oob=True).str.split("|")
        df = df.with_columns([
            first_ann.list.get(1, null_on_oob=True).alias("Effect"),
            first_ann.list.get(2, null_on_oob=True).alias("Impact"),
            first_ann.list.get(3, null_on_oob=True).alias("Gene"),
            first_ann.list.get(9, null_on_oob=True).alias("HGVS_c"),
            first_ann.list.get(10, null_on_oob=True).alias("HGVS_p"),
        ])
    else:
        # No annotation: add empty columns.
        # WHY: every report has the same columns, so scripts that read it don't break
        df = df.with_columns([
            pl.lit(None, dtype=pl.Utf8).alias(c)
            for c in ["Effect", "Impact", "Gene", "HGVS_c", "HGVS_p"]
        ])

    # Table 1: every call with its FILTER reason (the audit trail)
    all_calls = reports_dir / f"{sample}.all_calls.tsv"
    df.select(REPORT_COLS[:7] + ["FILTER"] + REPORT_COLS[7:]).write_csv(all_calls, separator="\t")

    # Table 2: only the calls that passed every filter
    passed = df.filter(pl.col("FILTER") == "PASS").select(REPORT_COLS)
    report = reports_dir / f"{sample}.high_confidence_variants.tsv"
    passed.write_csv(report, separator="\t")

    # Log how many calls each reason removed (a call can have several reasons)
    reasons = {r: df.filter(pl.col("FILTER").str.contains(r)).height for r in ["LowQual", "LowDepth", "LowAF"]}
    log.info(
        f"Filter QUAL >= {min_qual} & DP >= {min_dp} & AF >= {min_af}: "
        f"kept {passed.height} / {df.height} calls  (flagged: "
        + ", ".join(f"{k} {v}" for k, v in reasons.items()) + ")"
    )
    log.info(f"High-confidence report written → {report}")
    tmp_tsv.unlink(missing_ok=True)   # temporary file, no longer needed
    return report

# ---------------------------------------------------------------------------
# STEP 9b · STATISTICS – R (R/variant_stats.R)
#   WHAT: runs the R script on 06_reports/<sample>.all_calls.tsv and writes
#           06_reports/<sample>.r_stats.tsv     Ts/Tv with 95% CI, Wilcoxon PASS vs
#                                               flagged, Spearman QUAL~DP, AF spread,
#                                               KS test for even spread along the genome
#           06_reports/<sample>.r_windows.tsv   variants per window + Poisson/BH
#                                               hotspot test
#   WHY:  Polars (step 9) builds the tables; R is the natural tool for the
#         statistical tests. The dashboard (step 10) reads both outputs.
#   R runs from its own conda env `ngs-r`.
#   WHY a separate env: r-base needs a newer GSL library than bcftools 1.25,
#        so installing R into ngs-env downgrades bcftools to a broken 1.5.
# ---------------------------------------------------------------------------
R_SCRIPT       = Path(__file__).resolve().parent / "R" / "variant_stats.R"
HOTSPOT_WINDOW = 1000   # bp per window for the hotspot test


def find_rscript() -> Optional[str]:
    """Rscript to use: $NGS_RSCRIPT, then Rscript on PATH, then the ngs-r env."""
    import shutil
    envs = Path(sys.prefix).parent          # …/miniforge3/envs when run from ngs-env
    for c in (os.environ.get("NGS_RSCRIPT"), shutil.which("Rscript"),
              envs / "ngs-r" / "bin" / "Rscript",
              Path.home() / "miniforge3" / "envs" / "ngs-r" / "bin" / "Rscript"):
        if c and Path(c).exists():
            return str(c)
    return None


def run_r_stats(calls_tsv: Path, ref_fasta: Path, out_dir: Path, sample: str,
                window_bp: int = HOTSPOT_WINDOW) -> Tuple[Path, Path]:
    """Run R/variant_stats.R; returns (r_stats.tsv, r_windows.tsv)."""
    rscript = find_rscript()
    if rscript is None:
        raise RuntimeError("Rscript not found. Create the R env:  "
                           "mamba create -n ngs-r -c conda-forge r-base")
    log.info(f"Running R statistics ({Path(rscript).parent.parent.name})")
    # WHY the .fai: the R script needs the genome length for densities and windows
    fai = Path(str(ref_fasta) + ".fai")
    run_cmd([rscript, str(R_SCRIPT), str(calls_tsv), str(fai), str(out_dir), sample, str(window_bp)],
            out_dir / f"{sample}.r_stats.log")
    return out_dir / f"{sample}.r_stats.tsv", out_dir / f"{sample}.r_windows.tsv"

# ---------------------------------------------------------------------------
# MAIN – interactive run for one sample
#   WHAT: asks 4 questions, then runs steps 1–10 in order.
#   WHY:  the easiest way to run one sample; for many samples use Snakemake
#         (workflow/Snakefile), which calls the same step functions.
# ---------------------------------------------------------------------------
def main():
    print(f"\n{'='*65}")
    print(f"  NGS Pipeline {__version__}")
    print(f"  Author: {__author__}  <{__email__}>")
    print(f"{'='*65}\n")

    # ----- The 4 questions -----
    # WHY examples + re-asking: a typo here would otherwise only show up as a
    #     crash later, after downloads have already started.
    print("[1/4] Input reads – enter ONE of:")
    print("      • SRA run accession      e.g. SRR1553607")
    print("      • Local R1 FASTQ path    e.g. /workspaces/data/my_sample_R1.fastq.gz")
    print("                               (R2 is auto-detected: _R1→_R2 or _1→_2)")
    accession = ""
    while not accession:
        accession = input("      Accession or path: ").strip()
        if not accession:
            print("      ! This field is required.")

    print(f"\n[2/4] CPU threads – this machine has {os.cpu_count()} cores")
    print("      e.g. 4   (press Enter for default 8)")
    while True:
        t = input("      Threads [8]: ").strip()
        if not t:
            threads = 8
            break
        if t.isdigit() and int(t) > 0:
            threads = int(t)
            break
        print("      ! Please enter a positive whole number, e.g. 4")

    print("\n[3/4] Output folder – created if it doesn't exist")
    print("      e.g. ./NGS_run   or   ./SRR1553607_run   (press Enter for ./NGS_run)")
    base_dir = Path(input("      Output directory [./NGS_run]: ").strip() or "./NGS_run")

    # Optional for SRA (looked up from the organism), required for local files
    print("\n[4/4] Reference genome accession (RefSeq) – optional for SRA, REQUIRED for local FASTQ")
    print("      e.g. NC_002549.1 (Zaire ebolavirus)   NC_000913.3 (E. coli K-12)")
    print("           NC_045512.2 (SARS-CoV-2)         NC_003198.1 (Salmonella enterica)")
    print("      Press Enter to auto-select from the SRA organism name.")
    force_ref = input("      Reference accession: ").strip() or None
    print()

    # ----- STEP 1 · Setup -----
    dirs = setup_directories(base_dir)

    # ----- STEP 4 · Data input (local file or SRA download) -----
    if accession.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
        r1 = Path(accession)
        if not r1.exists():
            raise FileNotFoundError(f"Local FASTQ not found: {r1}")
        # Find R2 by changing only the read tag just before the extension:
        #   sample_R1.fastq.gz → sample_R2.fastq.gz   (also _1→_2, _R1_001→_R2_001)
        # WHY only at the end: changing every "_1" in the name turned
        #   sample_10_R1 into sample_20_R2, which could be another sample.
        mate_re = re.compile(r"(_R|_|R)1((?:_\d{3})?\.(?:fastq|fq)(?:\.gz)?)$")
        m = mate_re.search(r1.name)
        r2 = None
        if m:
            r2_path = r1.parent / mate_re.sub(r"\g<1>2\g<2>", r1.name)
            r2 = r2_path if r2_path.exists() else None
        organism = "unknown"   # local files carry no organism → reference must be given
        # Sample name = file name without _R1 and extension (my_sample_R1.fastq.gz → my_sample)
        # WHY: output files are then named after the sample
        sample = r1.name[:m.start()] if m else r1.name.split(".")[0]
        log.info(f"Using local FASTQ: {r1}" + (f" + {r2}" if r2 else " (single-end)"))
    else:
        organism, fastqs = process_accession(accession, dirs)
        sample = accession
        r1 = fastqs[0]
        r2 = fastqs[1] if len(fastqs) > 1 else None

    # ----- STEP 2 · Reference genome -----
    ref_fasta = download_reference(organism, dirs["ref"], force_accession=force_ref)

    # ----- STEP 3 · Index the reference -----
    index_reference(ref_fasta, dirs["logs"])

    # ----- STEP 5 · QC & trimming -----
    r1_trim, r2_trim = run_fastp(r1, r2, dirs["qc"], sample, threads)

    # ----- STEP 6 · Alignment -----
    bam = align_and_sort(ref_fasta, r1_trim, r2_trim, dirs["aligned"], sample, threads)

    # ----- STEP 7 · Variant calling -----
    raw_vcf = call_variants(ref_fasta, bam, dirs["variants"], sample)

    # ----- STEP 8 · Annotation -----
    # Run SnpEff only if SNPEFF_DB has a database for this reference.
    # WHY: a database for a different genome annotates nothing. Skipping it
    #      with a clear warning beats an empty result that looks real.
    # If SnpEff fails, continue with the raw VCF so you still get a report.
    snpeff_genome = SNPEFF_DB.get(ref_fasta.stem)
    if snpeff_genome is None:
        log.warning(
            f"No SnpEff database mapped for reference '{ref_fasta.stem}' – "
            "skipping annotation. To enable it, find the database with\n"
            f"    snpEff databases | grep -i <organism>\n"
            "and add  '<accession>': '<database>'  to SNPEFF_DB."
        )
        ann_vcf = raw_vcf
    else:
        try:
            ann_vcf = annotate_variants(raw_vcf, snpeff_genome, dirs["annotated"], sample)
        except Exception as e:
            log.warning(
                f"SnpEff annotation failed ({e}). Is the database downloaded?\n"
                f"    snpEff download -v {snpeff_genome}\n"
                "Using the raw VCF for reporting."
            )
            ann_vcf = raw_vcf

    # ----- STEP 9 · Filter & report -----
    report = filter_and_report(ann_vcf, dirs["reports"], sample)

    # ----- STEP 9b · Statistics in R -----
    # WHY only warn on failure: the report above is already complete; the R
    #     tables add statistics on top, and the dashboard works without them.
    try:
        run_r_stats(dirs["reports"] / f"{sample}.all_calls.tsv", ref_fasta, dirs["reports"], sample)
    except Exception as e:
        log.warning(f"R statistics skipped ({e}).")

    # ----- STEP 10 · Dashboard -----
    # One easy-to-read HTML page (charts + the CSV tables, viewable in the browser)
    # + the same CSV files for Power BI / Excel (see dashboard.py).
    # WHY: the TSV/VCF/BAM files are hard to read; the page shows the quality
    #      checks and the variants at a glance. If this step fails, the results
    #      above are still complete, so only warn.
    try:
        from dashboard import build_dashboard
        dash = build_dashboard(base_dir, min_qual=MIN_QUAL, min_dp=MIN_DP, min_af=MIN_AF)
    except Exception as e:
        log.warning(f"Dashboard not created ({e}). The report above is still complete.")
        dash = None

    print("\n" + "=" * 65)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print(f"  BAM           : {bam}")
    print(f"  Annotated VCF : {ann_vcf}")
    print(f"  Final report  : {report}")
    if dash:
        print(f"  Dashboard     : {dash}   (open in a browser)")
    print("=" * 65)

if __name__ == "__main__":
    main()