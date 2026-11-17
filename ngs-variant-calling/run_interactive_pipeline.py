#!/usr/bin/env python3
"""
NGS_Pipeline_v3.0 – Full automated variant calling pipeline
Author : Liam TrinhNguyen
Email  : LiamTrinhNguyen@gmail.com

Three ways to run it:
  1. Interactive, one sample:   python3 run_interactive_pipeline.py
  2. Snakemake, many samples:   snakemake -s workflow/Snakefile -c 4
  3. Nextflow, many samples:    bash workflow/run_nextflow.sh
     (2 and 3 both read their samples from config/config.yaml)

Steps 5b–5d (host removal, reference screen, Kraken2) run between QC and
alignment; see "Changes in v3.0" below.

The interactive script checks each answer as soon as it is typed: an SRA
accession is looked up on ENA (lookup_sra), a local FASTQ is summarised
(fastq_stats), and a reference accession is checked on NCBI
(check_reference). At question 4, "list" shows the built-in references and
"0" screens a read sample against all of them (run_reference_screen).

Changes in v3.0 (from v2.7). Sample checks before alignment, in all three
drivers (results: 02_qc/ and the dashboard's "Sample identity &
contamination" card). Variant calls are unchanged when the optional steps
are off (SRR1553607: 559 calls, 151 PASS – also with host removal on):
  - STEP 5b remove_host() (optional, host_index / NGS_HOST_INDEX): keeps only
    read pairs that don't map to a bwa-indexed host genome.
  - STEP 5c screen_sample() (on by default; screen: false / NGS_SCREEN=0):
    competitive mapping of 50,000 reads against every built-in genome plus
    human rRNA/mito decoys, reported per organism with shared reads split
    (one EM step), distant (85–95 %) reads, spread and EVENNESS (observed ÷
    expected spread). Lists every detected organism (mixed samples), flags
    divergent relatives, and warns when the configured reference isn't what
    the reads match. Option 0 at question 4 uses the same code.
  - STEP 5d classify_reads() (optional, kraken_db / NGS_KRAKEN_DB): Kraken2 +
    Bracken from the ngs-kraken env; falls back to Kraken2's species counts.
  - Negative controls (control: true in config.yaml): QC + checks only; the
    dashboard flags findings below 10× the controls' reads per million.
  - REF_SETS: multi-sequence references. Influenza A is now all 8 segments
    (Influenza_A_PR8); NC_002016.1 alone was segment 7 only.
  - sample_reads() retries a dropped ENA stream, and download_file() makes
    every FASTQ download resumable (HTTP Range, 5 attempts with back-off) and
    checks it against ENA's size and MD5; a verified file is reused. Snakemake
    (retries: 2) and Nextflow (errorStrategy 'retry') re-run failed downloads.

Changes in v2.7 (from v2.6). Interactive prompts only; the step functions and
their outputs are unchanged, so results match v2.6 (SRR1553607: 559 calls,
151 PASS):
  - Question 1 looks up an SRA accession on ENA straight away (organism,
    layout, reads, size, matching reference) and summarises a local FASTQ
    (reads, bases, lengths, GC, Q30). Typos, missing or broken files, and
    genome accessions (NC_…) are asked again instead of failing later.
  - Question 4 checks a typed accession on NCBI, rejects GCF_/GCA_
    assemblies, and adds "list" (numbered built-in references) and "0"
    (reference screen: competitive mapping of 50,000 reads, ranked by genome
    spread). It is required for local FASTQ.
  - REF_MAP moved to module level (reference_for() is shared with main()).
  - Piped answers that run out stop with a clear message (_ask()).

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
__version__ = "NGS_Pipeline_v3.0"

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

DOWNLOAD_TRIES = 5   # attempts per file (each resumes where the last one stopped)


def _md5(path: Path) -> str:
    import hashlib
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def download_file(url: str, dest: Path, size: Optional[int] = None, md5: Optional[str] = None,
                  tries: int = DOWNLOAD_TRIES) -> Path:
    """Download url to dest, RESUMING after a dropped connection and checking
    the size and MD5 that ENA publishes.
    WHY: ENA often closes long downloads mid-file ("Remote end closed connection",
         IncompleteRead). Starting again from zero can fail the same way for a
         big file; an HTTP Range request continues from the last byte received.
         A file that is already complete and correct is reused (fast re-runs)."""
    import time
    if dest.exists() and size and dest.stat().st_size == size and (not md5 or _md5(dest) == md5):
        log.info(f"{dest.name} already downloaded and verified")
        return dest
    last = None
    for attempt in range(1, tries + 1):
        have = dest.stat().st_size if dest.exists() else 0
        if size and have > size:                 # junk from an older, different file
            dest.unlink()
            have = 0
        try:
            if not (size and have == size):
                headers = {"Range": f"bytes={have}-"} if have else {}
                # timeout = (connect, between-bytes): a stalled stream fails instead of hanging
                with requests.get(url, stream=True, timeout=(30, 120), headers=headers) as r:
                    if have and r.status_code != 206:    # server ignored the range: start over
                        have = 0
                    r.raise_for_status()
                    total = size or (have + int(r.headers.get("content-length", 0)))
                    with open(dest, "ab" if have else "wb") as f, tqdm(
                        total=total, initial=have, unit="B", unit_scale=True, desc=dest.name, ncols=80
                    ) as bar:
                        for chunk in r.iter_content(chunk_size=1 << 20):
                            if chunk:
                                f.write(chunk)
                                bar.update(len(chunk))
            got = dest.stat().st_size
            if size and got != size:
                raise IOError(f"got {got:,} of {size:,} bytes")
            if md5 and _md5(dest) != md5:
                dest.unlink()                    # corrupt: the next attempt starts from zero
                raise IOError("MD5 checksum mismatch")
            return dest
        except Exception as e:                   # requests/urllib3 errors, short or corrupt file
            last = e
            if attempt < tries:
                wait = min(60, 5 * 2 ** (attempt - 1))
                log.warning(f"Download of {dest.name} interrupted ({type(e).__name__}: {str(e)[:80]}); "
                            f"resuming in {wait}s (attempt {attempt + 1}/{tries})")
                time.sleep(wait)
    raise RuntimeError(f"Could not download {dest.name} after {tries} attempts "
                       f"({type(last).__name__}). ENA may be busy – run again later; "
                       f"the partial file is kept and the next run resumes it.")


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
        "library_layout,library_strategy,read_count,base_count,fastq_ftp,fastq_bytes,fastq_md5"
    )
    log.info(f"Querying ENA for {run_accession}")
    import time
    for attempt in range(1, DOWNLOAD_TRIES + 1):
        try:
            resp = requests.get(url, timeout=90)
            if resp.status_code < 500:
                break
            err = f"HTTP {resp.status_code}"
        except requests.RequestException as e:
            err = type(e).__name__
        if attempt == DOWNLOAD_TRIES:
            raise RuntimeError(f"ENA API failed after {DOWNLOAD_TRIES} attempts ({err})")
        log.warning(f"ENA query failed ({err}); retrying in {5 * attempt}s")
        time.sleep(5 * attempt)
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
    links = data.get("fastq_ftp", "").split(";")
    sizes = data.get("fastq_bytes", "").split(";")
    md5s = data.get("fastq_md5", "").split(";")
    for i, link in enumerate(links):
        if not link.strip():
            continue
        clean = link.replace("ftp://", "").strip()
        file_url = f"https://{clean}"
        fname = out_dir / Path(file_url).name
        size = int(sizes[i]) if i < len(sizes) and sizes[i].isdigit() else None
        md5 = md5s[i].strip() if i < len(md5s) and md5s[i].strip() else None
        log.info(f"Downloading {fname.name}" + (f" ({size / 1e6:.1f} MB)" if size else "") + " …")
        download_file(file_url, fname, size=size, md5=md5)
        fastq_paths.append(fname)
        log.info(f"Saved {fname}" + (" (size and MD5 verified)" if md5 else ""))

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
# Organism name → RefSeq accession, used when no reference is given.
# Matching ignores case and works on part of the name ("Zaire ebolavirus"
# matches "Zaire ebolavirus Makona"); the first matching key wins.
# To support a new organism, add a line here (and to SNPEFF_DB for annotation).
# main() also uses this table: "list" at question 4 prints it.
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
    "Influenza A virus": "Influenza_A_PR8",     # all 8 segments, see REF_SETS

    # Eukaryotes (examples) – GCF_ assemblies can't be downloaded yet
    "Homo sapiens": "GCF_000001405.40",
    "Saccharomyces cerevisiae": "GCF_000146045.2",
}


# References made of several sequences: a name → the accessions downloaded
# together into ONE FASTA (<name>.fasta), used like any other reference
# (reference: Influenza_A_PR8 in config.yaml works too).
# WHY: influenza's genome is 8 separate RNA segments. The old single accession
#      NC_002016.1 was segment 7 only (1,027 of ~13,600 bp), so reads from the
#      other 7 segments were never aligned or called.
REF_SETS = {
    # Influenza A virus A/Puerto Rico/8/1934 (H1N1), segments 1–8
    "Influenza_A_PR8": ["NC_002023.1", "NC_002021.1", "NC_002022.1", "NC_002017.1",
                        "NC_002019.1", "NC_002018.1", "NC_002016.1", "NC_002020.1"],
}


def ref_accessions(ref: str) -> List[str]:
    """The NCBI accessions behind a reference name (one, or a REF_SETS entry)."""
    return REF_SETS.get(ref, [ref])


def reference_for(organism: str) -> Optional[Tuple[str, str]]:
    """(REF_MAP key, accession) for an organism name, or None if not mapped."""
    low = (organism or "").lower()
    key = next((k for k in REF_MAP if k.lower() in low), None)
    return (key, REF_MAP[key]) if key else None


def download_reference(organism: str, ref_dir: Path, force_accession: str = None) -> Path:
    """
    Download a suitable reference genome.
    You can force a specific accession with force_accession.
    """
    if force_accession:
        acc = force_accession.strip()
        log.info(f"Using forced reference accession: {acc}")
    else:
        match = reference_for(organism)
        key = match[0] if match else None

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
        # A REF_SETS name downloads all its sequences into one FASTA
        download_ncbi_fasta(",".join(ref_accessions(acc)), fasta)
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
# INPUT CHECKS – used by main() right after questions 1 and 4
#   WHAT: an SRA accession is looked up on ENA (organism, layout, size and the
#         reference it maps to); a local FASTQ is read once for its stats
#         (reads, bases, lengths, GC, quality); a reference is checked on NCBI.
#   WHY:  a typo or a wrong file shows up now, in seconds, instead of after
#         the downloads and alignment have already run.
# ---------------------------------------------------------------------------
FASTQ_EXT  = (".fastq", ".fastq.gz", ".fq", ".fq.gz")
SRA_PREFIX = ("SRR", "ERR", "DRR", "SAM", "PRN", "SRS", "ERS", "DRS")
# R2 = R1 with only the read tag before the extension changed (_R1→_R2, _1→_2, R1_001→R2_001)
MATE_RE    = re.compile(r"(_R|_|R)1((?:_\d{3})?\.(?:fastq|fq)(?:\.gz)?)$")
# Files bigger than this are sampled, not read to the end (a 300 MB .gz takes ~20 s)
SCAN_LIMIT_BYTES  = 300 * 1024 ** 2
SCAN_SAMPLE_READS = 1_000_000


def _human_bytes(n: float) -> str:
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024


def lookup_sra(accession: str, timeout: int = 30) -> dict:
    """Ask ENA about a run without downloading anything.
    Raises LookupError if ENA has no such run, ConnectionError if ENA can't be reached."""
    if requests is None:
        raise ConnectionError("the 'requests' library is not installed")
    fields = ("run_accession,scientific_name,library_layout,library_strategy,"
              "instrument_platform,instrument_model,read_count,base_count,"
              "fastq_bytes,fastq_ftp,study_accession")
    try:
        r = requests.get("https://www.ebi.ac.uk/ena/portal/api/filereport",
                         params={"accession": accession, "result": "read_run", "fields": fields},
                         timeout=timeout)
    except requests.RequestException as e:
        raise ConnectionError(type(e).__name__) from e
    if r.status_code == 400:            # ENA's answer to a malformed accession
        raise LookupError(f"ENA doesn't recognise '{accession}' as a read accession")
    if r.status_code != 200:
        raise ConnectionError(f"ENA answered HTTP {r.status_code}")
    lines = r.text.strip().splitlines()
    if len(lines) < 2:                  # header only: well-formed but unknown
        raise LookupError(f"ENA has no public run called '{accession}'")
    header = lines[0].split("\t")
    info = dict(zip(header, lines[1].split("\t")))
    sizes = [int(b) for b in info.get("fastq_bytes", "").split(";") if b.isdigit()]
    info.update(n_runs=len(lines) - 1, n_files=len(sizes), total_bytes=sum(sizes))
    return info


def print_sra_info(accession: str, info: dict) -> None:
    """The ENA answer as a short card under question 1."""
    organism = info.get("scientific_name") or "unknown"
    layout   = info.get("library_layout") or "?"
    reads, bases = info.get("read_count", ""), info.get("base_count", "")
    model    = info.get("instrument_model") or info.get("instrument_platform") or "?"
    match    = reference_for(organism)
    if match and not match[1].startswith("GCF_"):
        ref = f"{match[1]} ({match[0]}) – used if you press Enter at question 4"
    elif match:
        ref = f"{match[1]} – an assembly that can't be downloaded yet; pick one at question 4"
    else:
        ref = "none built in – give one at question 4 (type  list  there to see them)"
    rows = [
        ("Organism", organism),
        ("Layout", layout + {"PAIRED": " (R1 + R2)", "SINGLE": " (one read per fragment)"}.get(layout, "")),
        ("Reads", (f"{int(reads):,}" if reads.isdigit() else "?") + f" · {model}"),
        ("Bases", f"{int(bases):,} bp" if bases.isdigit() else "?"),
        ("Library", info.get("library_strategy") or "?"),
        ("Download", f"{_human_bytes(info['total_bytes'])} in {info['n_files']} file(s)" if info["n_files"] else "no FASTQ files listed"),
        ("Study", info.get("study_accession") or "?"),
        ("Reference", ref),
    ]
    print(f"\n      ✓ Found on ENA: {info.get('run_accession') or accession}")
    for k, v in rows:
        print(f"        {k:<10} {v}")
    if info["n_runs"] > 1:
        print(f"        ! {accession} covers {info['n_runs']} runs; only the first ({info.get('run_accession')}) is used.")
    if layout == "SINGLE":
        print("        Note: single-end runs work here, but Snakemake/Nextflow need paired-end reads.")
    if not info["n_files"]:
        print("        ! ENA lists no FASTQ files for this run, so the download in step 4 will fail.")
    print()


def fastq_stats(path: Path) -> dict:
    """Read a FASTQ (.gz or plain) once: reads, bases, read lengths, GC and quality.
    Files over SCAN_LIMIT_BYTES are sampled (first SCAN_SAMPLE_READS reads) and
    the totals estimated from the share of the file that was read.
    Raises ValueError if the file isn't valid FASTQ."""
    import gzip
    size = path.stat().st_size
    limit = SCAN_SAMPLE_READS if size > SCAN_LIMIT_BYTES else None
    low_q = bytes(range(33, 63))        # Phred+33 characters for Q0–Q29
    n = bases = gc = q30 = qsum = 0
    lmin, lmax, first_id = None, 0, ""
    with open(path, "rb") as raw:
        f = gzip.GzipFile(fileobj=raw) if path.name.endswith(".gz") else raw
        try:
            while limit is None or n < limit:
                head = f.readline()
                if not head:
                    break
                seq, plus, qual = f.readline().rstrip(), f.readline(), f.readline().rstrip()
                if not head.startswith(b"@") or not plus.startswith(b"+") or len(qual) != len(seq):
                    raise ValueError(f"record {n + 1} isn't valid FASTQ (line {4 * n + 1})")
                if n == 0:
                    first_id = head[1:].split()[0].decode(errors="replace") if head[1:].split() else ""
                L = len(seq)
                n += 1
                bases += L
                lmin = L if lmin is None else min(lmin, L)
                lmax = max(lmax, L)
                gc += seq.count(b"G") + seq.count(b"C") + seq.count(b"g") + seq.count(b"c")
                q30 += len(qual.translate(None, low_q))
                qsum += sum(qual) - 33 * L
            sampled = limit is not None and n >= limit and bool(f.readline())
        except (OSError, EOFError) as e:     # gzip: not a .gz file, or cut short
            raise ValueError(f"can't read it as gzip ({e})") from e
        share = raw.tell() / size if sampled and size else 1.0
    if n == 0:
        raise ValueError("the file has no reads")
    return {
        "size": size, "reads": n, "bases": bases, "sampled": sampled,
        "est_reads": int(n / share), "est_bases": int(bases / share),
        "min_len": lmin, "max_len": lmax, "mean_len": bases / n,
        "gc_pct": 100 * gc / bases if bases else 0.0,
        "q30_pct": 100 * q30 / bases if bases else 0.0,
        "mean_q": qsum / bases if bases else 0.0,
        "first_id": first_id,
    }


def print_fastq_info(label: str, path: Path, st: dict) -> None:
    """The stats of one FASTQ file as a short card under question 1."""
    print(f"      ✓ {label}: {path.name} ({_human_bytes(st['size'])})")
    if st["sampled"]:
        print(f"        Reads      ≈ {st['est_reads']:,}   (estimated from the first {st['reads']:,} reads)")
        print(f"        Bases      ≈ {st['est_bases']:,} bp")
    else:
        print(f"        Reads      {st['reads']:,}")
        print(f"        Bases      {st['bases']:,} bp")
    lengths = f"{st['min_len']}" if st["min_len"] == st["max_len"] else f"{st['min_len']}–{st['max_len']}"
    print(f"        Length     {lengths} bp (mean {st['mean_len']:.1f})")
    print(f"        GC         {st['gc_pct']:.1f} %")
    print(f"        Quality    mean Q{st['mean_q']:.1f} · {st['q30_pct']:.1f} % of bases ≥ Q30")
    if st["first_id"]:
        print(f"        First read {st['first_id']}")


def describe_local_fastq(r1: Path) -> bool:
    """Stats for R1 (and R2, if found next to it). False if R1 isn't valid FASTQ."""
    m = MATE_RE.search(r1.name)
    r2 = r1.parent / MATE_RE.sub(r"\g<1>2\g<2>", r1.name) if m else None
    sample = r1.name[:m.start()] if m else r1.name.split(".")[0]
    print("      Reading the file for a summary …", flush=True)
    try:
        st1 = fastq_stats(r1)
    except ValueError as e:
        print(f"      ! {r1.name}: {e}.")
        return False
    print()
    print_fastq_info("R1", r1, st1)
    if r2 is not None and r2.is_file():
        try:
            st2 = fastq_stats(r2)
            print_fastq_info("R2", r2, st2)
            if not (st1["sampled"] or st2["sampled"]) and st1["reads"] != st2["reads"]:
                print(f"        ! R1 has {st1['reads']:,} reads but R2 has {st2['reads']:,}: "
                      "they may not be a pair.")
        except ValueError as e:
            print(f"      ! R2 {r2.name}: {e}. Fix or remove it before running.")
    elif r2 is not None:
        print(f"      No R2 found (looked for {r2.name}): runs as single-end.")
    else:
        print("      No _R1/_1 tag before the extension, so no R2 is looked for: runs as single-end.")
    print(f"      Sample name: {sample}  (used in the output file names)")
    print("      Local files carry no organism name, so question 4 needs a reference.\n")
    return True


def check_reference(accession: str) -> dict:
    """Title and total length of a reference on NCBI (a REF_SETS name checks
    every sequence in it). Raises LookupError if NCBI has no such record,
    ConnectionError if it can't be reached."""
    if Entrez is None:
        raise ConnectionError("Biopython is not installed")
    ids = ref_accessions(accession)
    try:
        handle = Entrez.esummary(db="nuccore", id=",".join(ids))
        record = Entrez.read(handle)
        handle.close()
    except RuntimeError as e:            # Entrez: "Invalid uid …"
        raise LookupError(f"NCBI Nucleotide has no record '{accession}'") from e
    except Exception as e:               # network errors (URLError, timeouts …)
        raise ConnectionError(type(e).__name__) from e
    if len(record) != len(ids):
        raise LookupError(f"NCBI returned {len(record)} of the {len(ids)} sequences of '{accession}'")
    length = sum(int(d.get("Length", 0)) for d in record)
    title = str(record[0].get("Title", ""))
    if len(ids) > 1:
        title = f"{len(ids)} sequences, e.g. {title}"
    return {"title": title, "length": length}


def reference_choices() -> List[Tuple[str, str]]:
    """Each built-in reference once: (accession, organism name)."""
    seen = {}
    for name, acc in REF_MAP.items():
        seen.setdefault(acc, name)
    return list(seen.items())


def print_reference_list() -> None:
    """The numbered table shown when you type  list  at question 4."""
    print("\n      Built-in references – type the number or the accession:")
    print("         0. all                Screen your reads against every reference below and")
    print("                               suggest the best match (downloads them once, ~17 MB)")
    for i, (acc, name) in enumerate(reference_choices(), 1):
        if acc in SNPEFF_DB:
            note = "with gene annotation (SnpEff)"
        elif acc in REF_SETS:
            note = f"{len(REF_SETS[acc])} sequences (all segments)"
        elif acc.startswith("GCF_"):
            note = "can't be downloaded yet (assembly)"
        else:
            note = ""
        print(f"      {i:>4}. {acc:<18} {name:<28} {note}")
    print("      Not listed? Type any RefSeq nucleotide accession (search https://www.ncbi.nlm.nih.gov/nuccore).\n")


# ---------------------------------------------------------------------------
# STEP 5b · HOST REMOVAL (optional) – drop reads from the host (e.g. human)
#   WHAT: maps the trimmed reads to a host bwa index and keeps only pairs in
#         which NEITHER read maps → 02_qc/<sample>_R1/_R2.host_removed.fastq.gz
#         + <sample>.host_removal.json (reads in, kept, host %). Later steps
#         use the kept reads.
#   WHY:  clinical samples are often > 90 % host. Removing those reads speeds
#         up every later step, stops host reads aligning to conserved pathogen
#         genes (false variants), and keeps human sequence out of the outputs.
#   HOW TO TURN IT ON: a bwa index of the host genome –
#         host_index: /path/T2T-CHM13.fa   in config.yaml (Snakemake/Nextflow)
#         export NGS_HOST_INDEX=/path/T2T-CHM13.fa   for the interactive script
#         (build once with  bwa index T2T-CHM13.fa : ~1 h and ~5 GB RAM).
#         Empty (default) = step skipped.
# ---------------------------------------------------------------------------
HOST_INDEX = os.environ.get("NGS_HOST_INDEX", "")


def remove_host(r1: Path, r2: Optional[Path], host_index: str, out_dir: Path,
                sample: str, threads: int) -> Tuple[Path, Optional[Path]]:
    """Keep only reads (pairs) that don't map to the host. Returns the new R1/R2."""
    import shlex
    idx = Path(host_index).expanduser()        # "~/host/chm13.fa" from a config file
    if not Path(str(idx) + ".bwt").exists():
        raise RuntimeError(f"No bwa index at {idx} (expected {idx}.bwt). Build it with:  bwa index {idx}")
    bam = out_dir / f"{sample}.host_check.bam"
    o1 = out_dir / f"{sample}_R1.host_removed.fastq.gz"
    o2 = out_dir / f"{sample}_R2.host_removed.fastq.gz" if r2 else None
    reads = " ".join(shlex.quote(str(p)) for p in ([r1, r2] if r2 else [r1]))
    log.info(f"Removing host reads ({idx.name})")
    run_cmd(["bash", "-c", f"set -o pipefail; bwa mem -t {threads} {shlex.quote(str(idx))} {reads} "
                           f"2> {shlex.quote(str(out_dir / f'{sample}.host_removal.log'))} "
                           f"| samtools view -b -o {shlex.quote(str(bam))} -"])
    # Paired: keep a pair only if BOTH reads are unmapped (-f 12); single-end: -f 4
    keep = "12" if r2 else "4"
    total = int(run_cmd(["samtools", "view", "-c", "-F", "0x900", str(bam)]).stdout.strip())
    kept = int(run_cmd(["samtools", "view", "-c", "-f", keep, "-F", "0x900", str(bam)]).stdout.strip())
    if r2:
        run_cmd(["samtools", "fastq", "-f", keep, "-F", "0x900", "-n", "-1", str(o1), "-2", str(o2),
                 "-0", "/dev/null", "-s", "/dev/null", str(bam)])
    else:
        run_cmd(["samtools", "fastq", "-f", keep, "-F", "0x900", "-n", "-0", str(o1), str(bam)])
    bam.unlink()
    stats = {"host_index": str(idx), "reads_in": total, "reads_kept": kept,
             "host_reads": total - kept, "host_pct": round(100 * (total - kept) / total, 2) if total else 0.0}
    (out_dir / f"{sample}.host_removal.json").write_text(json.dumps(stats, indent=2))
    log.info(f"Host reads removed: {stats['host_reads']:,} of {total:,} ({stats['host_pct']} %)")
    return o1, o2


# ---------------------------------------------------------------------------
# STEP 5d · CLASSIFICATION (optional) – what is in the sample? Kraken2 + Bracken
#   WHAT: Kraken2 assigns every read to a taxon using a large database (an
#         OPEN-world answer, unlike the reference screen's fixed list); Bracken
#         re-estimates the number of reads per species.
#         → 02_qc/<sample>.kraken2.report, <sample>.bracken_species.tsv
#   WHY:  the right first step for an unknown sample, and the way to see
#         contamination. Compare with a negative control (control: true).
#   HOW TO TURN IT ON: a Kraken2 database folder (hash.k2d, …) –
#         kraken_db: /path/k2_standard_08gb   in config.yaml
#         export NGS_KRAKEN_DB=/path/k2_standard_08gb   for the interactive script
#         Prebuilt databases: https://benlangmead.github.io/aws-indexes/k2
#         (Standard-8: 8 GB RAM; Viral: 1 GB, viruses only). A database larger
#         than RAM needs --memory-mapping (kraken_memory_mapping: true /
#         NGS_KRAKEN_MMAP=1): slower, but it runs.
#         Kraken2/Bracken live in their own env, ngs-kraken (like R in ngs-r).
# ---------------------------------------------------------------------------
KRAKEN_DB = os.environ.get("NGS_KRAKEN_DB", "")
KRAKEN_MMAP = os.environ.get("NGS_KRAKEN_MMAP", "") not in ("", "0", "false", "False")
BRACKEN_COLUMNS = ["name", "taxonomy_id", "taxonomy_lvl", "kraken_assigned_reads",
                   "added_reads", "new_est_reads", "fraction_total_reads"]


def find_env_tool(name: str, env: str) -> Optional[str]:
    """A tool from PATH, else from a sibling conda env (…/envs/<env>/bin/<name>)."""
    import shutil
    envs = Path(sys.prefix).parent
    for c in (shutil.which(name), envs / env / "bin" / name,
              Path.home() / "miniforge3" / "envs" / env / "bin" / name):
        if c and Path(c).exists():
            return str(c)
    return None


def _mean_read_length(fastq: Path, n: int = 2000) -> float:
    import gzip
    opener = gzip.open if fastq.name.endswith(".gz") else open
    total = count = 0
    with opener(fastq, "rb") as fh:
        for i, line in enumerate(fh):
            if i % 4 == 1:
                total += len(line.rstrip())
                count += 1
                if count >= n:
                    break
    return total / count if count else 0.0


def _species_from_report(report: Path, out: Path) -> None:
    """Fallback when Bracken can't run: species rows (rank S) of the Kraken2
    report, written with Bracken's columns (new_est_reads = Kraken's clade reads)."""
    rows = []
    for line in report.read_text().splitlines():
        f = line.split("\t")
        if len(f) >= 6 and f[3] == "S":
            rows.append([f[5].strip(), f[4], "S", int(f[1]), 0, int(f[1])])
    total = sum(r[5] for r in rows) or 1
    with open(out, "w") as fh:
        fh.write("\t".join(BRACKEN_COLUMNS) + "\n")
        for r in sorted(rows, key=lambda r: -r[5]):
            fh.write("\t".join(map(str, r + [f"{r[5] / total:.5f}"])) + "\n")


def classify_reads(r1: Path, r2: Optional[Path], db: str, out_dir: Path, sample: str,
                   threads: int, memory_mapping: bool = KRAKEN_MMAP) -> Tuple[Path, Path]:
    """Kraken2 (+ Bracken) on the reads. Returns (kraken2 report, species table)."""
    db = str(Path(db).expanduser())            # "~/kraken_db/viral" from a config file
    kraken = find_env_tool("kraken2", "ngs-kraken")
    if kraken is None:
        raise RuntimeError("kraken2 not found. Create its env:  "
                           "mamba create -n ngs-kraken -c conda-forge -c bioconda -y kraken2 bracken")
    if not (Path(db) / "hash.k2d").exists():
        raise RuntimeError(f"No Kraken2 database at {db} (hash.k2d missing)")
    report = out_dir / f"{sample}.kraken2.report"
    species = out_dir / f"{sample}.bracken_species.tsv"
    cmd = [kraken, "--db", str(db), "--threads", str(threads),
           "--report", str(report), "--output", "/dev/null"]
    if memory_mapping:
        cmd.append("--memory-mapping")
    cmd += ["--paired", str(r1), str(r2)] if r2 else [str(r1)]
    log.info(f"Classifying reads with Kraken2 ({Path(db).name})")
    run_cmd(cmd, out_dir / f"{sample}.kraken2.log")

    # Bracken needs a k-mer distribution built for (about) this read length
    bracken = find_env_tool("bracken", "ngs-kraken")
    lengths = sorted(int(m.group(1)) for f in Path(db).glob("database*mers.kmer_distrib")
                     if (m := re.match(r"database(\d+)mers", f.name)))
    read_len = _mean_read_length(r1)
    use = max((L for L in lengths if L <= read_len + 5), default=None)
    try:
        if bracken is None or use is None:
            raise RuntimeError("Bracken or a matching kmer_distrib file is missing")
        run_cmd([bracken, "-d", str(db), "-i", str(report), "-o", str(species),
                 "-w", str(out_dir / f"{sample}.bracken.report"), "-r", str(use), "-l", "S", "-t", "10"],
                out_dir / f"{sample}.bracken.log")
        if not species.exists():
            raise RuntimeError("Bracken wrote no output")
    except RuntimeError as e:
        log.warning(f"Bracken skipped ({e}); using Kraken2's species counts instead")
        _species_from_report(report, species)
    return report, species


# ---------------------------------------------------------------------------
# STEP 5c · REFERENCE SCREEN – which genome do the reads come from?
#   WHAT: maps a sample of the reads (SCREEN_READS, 50,000) with bwa mem against
#         ALL downloadable built-in references at once, plus human rRNA and
#         mitochondrial DNA as host decoys, and reports per ORGANISM (all
#         sequences of a genome, e.g. 8 influenza segments, count together):
#           Reads      reads placed on this organism only (strict: ≥ 80 % of the
#                      read aligned, ≤ 5 % mismatches)
#           Shared     reads that fit several organisms equally well, split in
#                      proportion to each one's own reads (one EM step) instead
#                      of being thrown away
#           Distant    reads aligning at 85–95 % identity: a divergent strain or
#                      a relative that is not in the list
#           Spread     how many of SCREEN_WINDOWS (100) equal stretches along the
#                      genome got a uniquely placed (MAPQ ≥ 1) strict read
#           Evenness   Spread ÷ the spread expected if those reads landed at random,
#                      100·(1 − e^(−reads/100)). ≈ 1: reads all along the genome
#                      (the organism is there); ≪ 1: piled up in a few genes
#           Status     detected (Spread ≥ 5 and Evenness ≥ 0.5) / divergent? (the
#                      same test on the 85–95 % reads) / piled up / host / –
#   Used two ways: option "0" at question 4 (pick a reference), and as a
#   pipeline step for every sample (screen_sample → 02_qc/<sample>.reference_screen.tsv)
#   that warns when the configured reference is not what the reads match.
#   WHY spread and evenness, not counts: on SRR1553607 (Ebola) 68 strict reads
#        hit M. tuberculosis, all inside its rRNA genes (conserved bacterial rRNA,
#        e.g. from a reagent contaminant), versus 39 Ebola reads spread over
#        32/100 of the virus (32 expected: evenness 1.0). Spread alone was not
#        enough: S. aureus carries several rRNA copies around its genome (6/100
#        windows from 59 reads, evenness 0.13), and conserved genes of related
#        contaminants gave M. tuberculosis 211 reads at 85–95 % over 15/100
#        windows (evenness 0.17). Every detected organism is listed, so mixed
#        infections show.
#   LIMITS: a closed list – it can only rank these genomes. For truly unknown
#        samples use classify_reads() (Kraken2) and a negative control.
#   The panel is downloaded and indexed once into SCREEN_CACHE (a file lock
#   stops parallel Snakemake/Nextflow jobs from building it twice).
# ---------------------------------------------------------------------------
SCREEN_READS = 50_000
SCREEN_CACHE = Path(os.environ.get("NGS_SCREEN_CACHE",
                                   Path.home() / ".cache" / "ngs_pipeline" / "reference_screen"))
SCREEN_WINDOWS    = 100      # each genome is split into this many equal stretches
SCREEN_MIN_SPREAD = 5        # an organism needs reads in at least this many to count as detected
SCREEN_PILED_READS = 20      # this many reads in ≤ SCREEN_MIN_SPREAD stretches = "piled up"
# Host decoys: give human reads a better home than a conserved pathogen gene
SCREEN_HOST = ("Human (rRNA + mitochondria)", ["U13369.1", "NC_012920.1"])
SCREEN_MIN_EVENNESS = 0.5    # observed ÷ expected spread needed to count as detected
SCREEN_COLUMNS = ["Organism", "Reference", "Kind", "Reads", "Shared_reads", "Distant_reads",
                  "Spread", "Evenness", "Distant_spread", "Distant_evenness", "Mean_identity",
                  "Reads_per_million", "Sampled_reads", "Status"]


def _expected_spread(placed_reads: int) -> float:
    """Windows expected to get ≥ 1 read if placed_reads landed at random."""
    import math
    return SCREEN_WINDOWS * (1 - math.exp(-placed_reads / SCREEN_WINDOWS))


def screen_panel() -> List[Tuple[str, str, str, List[str]]]:
    """(organism, reference name, kind, accessions) for every screened genome."""
    panel = [(name, ref, "target", ref_accessions(ref))
             for ref, name in reference_choices() if not ref.startswith(("GCF_", "GCA_"))]
    panel.append((SCREEN_HOST[0], "host", "host", SCREEN_HOST[1]))
    return panel


def screen_index() -> Path:
    """Combined FASTA of the whole panel, bwa-indexed, cached in SCREEN_CACHE."""
    import fcntl
    panel = screen_panel()
    wanted = [f"{acc}\t{ref}" for _, ref, _, accs in panel for acc in accs]
    SCREEN_CACHE.mkdir(parents=True, exist_ok=True)
    fasta, listed = SCREEN_CACHE / "all_references.fasta", SCREEN_CACHE / "panel.tsv"
    with open(SCREEN_CACHE / ".lock", "w") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)          # one builder at a time
        if (fasta.exists() and listed.exists() and Path(str(fasta) + ".bwt").exists()
                and listed.read_text().splitlines() == wanted):
            return fasta
        if Entrez is None:
            raise RuntimeError("Biopython is needed to download the references (pip install biopython)")
        accs = [line.split("\t")[0] for line in wanted]
        log.info(f"Reference screen: downloading {len(accs)} sequences from NCBI (first time only)")
        handle = Entrez.efetch(db="nuccore", id=",".join(accs), rettype="fasta", retmode="text")
        data = handle.read()
        handle.close()
        got = {line[1:].split()[0] for line in data.splitlines() if line.startswith(">")}
        missing = [a for a in accs if a not in got]
        if missing:
            raise RuntimeError(f"NCBI didn't return {', '.join(missing)}")
        fasta.write_text(data)
        log.info("Reference screen: indexing with bwa (about a minute, first time only)")
        run_cmd(["bwa", "index", str(fasta)], SCREEN_CACHE / "bwa_index.log")
        listed.write_text("\n".join(wanted) + "\n")
    return fasta


def sample_reads(source, out_fq: Path, n: int = SCREEN_READS, tries: int = 3) -> int:
    """Copy the first n reads of a FASTQ into out_fq. source is a local Path or
    an https URL of a .fastq.gz (streamed: only the start is downloaded).
    A dropped connection is retried (tries times), then raised as ConnectionError."""
    import gzip
    import time

    def copy(stream) -> int:
        count = 0
        with open(out_fq, "wb") as out:
            while count < n:
                rec = [stream.readline() for _ in range(4)]
                if not rec[3]:
                    break
                out.writelines(rec)
                count += 1
        return count

    if isinstance(source, Path):
        with open(source, "rb") as raw:
            return copy(gzip.GzipFile(fileobj=raw) if source.name.endswith(".gz") else raw)
    if requests is None:
        raise ConnectionError("the 'requests' library is not installed")
    last = None
    for attempt in range(1, tries + 1):
        try:
            with requests.get(source, stream=True, timeout=60) as resp:
                resp.raise_for_status()
                resp.raw.decode_content = True
                return copy(gzip.GzipFile(fileobj=resp.raw) if source.endswith(".gz") else resp.raw)
        except Exception as e:      # requests/urllib3 errors, a cut-off gzip stream …
            last = e
            log.warning(f"Streaming reads failed ({type(e).__name__}), attempt {attempt}/{tries}")
            time.sleep(2 * attempt)
    raise ConnectionError(type(last).__name__)


def _aligned_len(cigar: str) -> int:
    return sum(int(n) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", cigar) if op in "M=X")


def screen_reads(reads_fq: Path, threads: int) -> Tuple[List[dict], int]:
    """Map a read sample against the panel; one result row per organism
    (columns as SCREEN_COLUMNS), best first. Returns (rows, reads sampled)."""
    index_fa = screen_index()
    panel = screen_panel()
    group_of = {acc: ref for _, ref, _, accs in panel for acc in accs}
    info = {ref: (name, kind) for name, ref, kind, _ in panel}
    lengths, offsets, glen = {}, {}, {}
    unique, distant, shared_sets, ident = {}, {}, [], {}
    placed, placed_d = {}, {}                 # MAPQ ≥ 1 reads: the ones that count for spread
    win_close, win_dist = {}, {}
    n_reads = sum(1 for _ in open(reads_fq)) // 4

    def window(contig, pos):
        g = group_of[contig]
        return g, (offsets[contig] + pos - 1) * SCREEN_WINDOWS // glen[g]

    proc = subprocess.Popen(["bwa", "mem", "-t", str(threads), str(index_fa), str(reads_fq)],
                            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True)
    for line in proc.stdout:
        if line.startswith("@"):
            if line.startswith("@SQ"):
                tags = dict(t.split(":", 1) for t in line.rstrip("\n").split("\t")[1:])
                lengths[tags["SN"]] = int(tags["LN"])
            continue
        if not offsets:                       # first record: sequence lengths are known
            for contig, ln in lengths.items():
                g = group_of.get(contig, contig)
                group_of.setdefault(contig, g)
                offsets[contig] = glen.get(g, 0)
                glen[g] = glen.get(g, 0) + ln
        f = line.rstrip("\n").split("\t")
        flag, contig, pos, mapq, cigar = int(f[1]), f[2], int(f[3]), int(f[4]), f[5]
        if flag & 0x904:                      # unmapped, secondary or supplementary
            continue
        aligned = _aligned_len(cigar)
        if aligned < 0.8 * len(f[9]):         # mostly unaligned: a short hit on a conserved motif
            continue
        nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), 0)
        identity = 1 - nm / aligned
        g, w = window(contig, pos)
        # Other places the read fits as well (bwa's XA tag) → which organisms are candidates?
        cands = {g}
        if mapq == 0:
            xa = next((t[5:] for t in f[11:] if t.startswith("XA:Z:")), "")
            for hit in filter(None, xa.split(";")):
                c, _, cig, nm2 = hit.split(",")
                if 1 - int(nm2) / max(1, _aligned_len(cig)) >= identity - 0.01:
                    cands.add(group_of.get(c, c))
        if identity >= 0.95:
            if len(cands) == 1:               # one organism (several copies inside it are fine)
                unique[g] = unique.get(g, 0) + 1
                ident[g] = ident.get(g, 0.0) + identity
                if mapq >= 1:                 # WHY: a repeat copy's position is a guess
                    placed[g] = placed.get(g, 0) + 1
                    win_close.setdefault(g, set()).add(w)
            else:
                shared_sets.append(cands)
        elif identity >= 0.85 and len(cands) == 1:
            distant[g] = distant.get(g, 0) + 1
            if mapq >= 1:
                placed_d[g] = placed_d.get(g, 0) + 1
                win_dist.setdefault(g, set()).add(w)
    if proc.wait() != 0:
        raise RuntimeError("bwa mem failed during the reference screen")

    # One EM step: split each shared read between its candidates in proportion
    # to their own (unique) reads; evenly if none of them has any
    shared = {}
    for cands in shared_sets:
        weights = {c: unique.get(c, 0) for c in cands}
        total = sum(weights.values())
        for c in cands:
            shared[c] = shared.get(c, 0.0) + (weights[c] / total if total else 1 / len(cands))

    rows = []
    for name, ref, kind, _ in panel:
        u, d = unique.get(ref, 0), distant.get(ref, 0)
        spread, dspread = len(win_close.get(ref, ())), len(win_dist.get(ref, ()))
        even = spread / _expected_spread(placed[ref]) if placed.get(ref) else 0.0
        deven = dspread / _expected_spread(placed_d[ref]) if placed_d.get(ref) else 0.0
        if kind == "host":
            status = "host" if u + shared.get(ref, 0) > 0 else "–"
        elif spread >= SCREEN_MIN_SPREAD and even >= SCREEN_MIN_EVENNESS:
            status = "detected"
        elif dspread >= SCREEN_MIN_SPREAD and deven >= SCREEN_MIN_EVENNESS:
            status = "divergent?"
        elif max(u, d) >= SCREEN_PILED_READS:
            status = "piled up"
        else:
            status = "–"
        reads = u + shared.get(ref, 0)
        rows.append({
            "Organism": name, "Reference": ref, "Kind": kind, "Reads": u,
            "Shared_reads": round(shared.get(ref, 0.0), 1), "Distant_reads": d,
            "Spread": spread, "Evenness": round(even, 2),
            "Distant_spread": dspread, "Distant_evenness": round(deven, 2),
            "Mean_identity": round(100 * ident[ref] / u, 2) if u else None,
            "Reads_per_million": round(1e6 * reads / n_reads, 1) if n_reads else 0.0,
            "Sampled_reads": n_reads, "Status": status,
        })
    order = {"detected": 0, "divergent?": 1, "piled up": 2, "host": 3, "–": 4}
    rows.sort(key=lambda r: (order[r["Status"]], -r["Spread"], -r["Distant_spread"], -r["Reads"]))
    return rows, n_reads


def screen_findings(rows: List[dict], expected_ref: Optional[str] = None) -> List[str]:
    """Plain-language warnings about a screen result (also logged by screen_sample)."""
    notes = []
    detected = [r for r in rows if r["Status"] == "detected"]
    n = rows[0]["Sampled_reads"] if rows else 0
    if not detected:
        notes.append("No built-in genome is covered well: the organism is probably not in the "
                     "list (classify the reads with Kraken2, or give its accession).")
    if len(detected) > 1:
        notes.append("More than one organism detected (" +
                     ", ".join(f"{r['Organism']} {r['Spread']}/{SCREEN_WINDOWS}" for r in detected) +
                     "): possible mixed infection or contamination – compare with a negative control.")
    for r in rows:
        if r["Status"] == "divergent?":
            notes.append(f"{r['Organism']}: {r['Distant_reads']} reads at 85–95 % identity spread over "
                         f"{r['Distant_spread']}/{SCREEN_WINDOWS} – a divergent strain or a relative not in "
                         "the list. Variant calls against this reference would be unreliable; "
                         "consider a closer reference or de novo assembly.")
    host = next((r for r in rows if r["Kind"] == "host"), None)
    if host and n and host["Reads"] + host["Shared_reads"] > 0.2 * n:
        notes.append(f"{100 * (host['Reads'] + host['Shared_reads']) / n:.0f} % of the sampled reads are "
                     "human rRNA/mitochondrial: consider host removal (host_index).")
    if expected_ref:
        exp = next((r for r in rows if r["Reference"] == expected_ref), None)
        if exp is None:
            pass                              # reference not in the panel: nothing to compare
        elif exp["Status"] != "detected":
            notes.append(f"The configured reference {expected_ref} is NOT detected in the reads "
                         f"(spread {exp['Spread']}/{SCREEN_WINDOWS}). Check the sample and the reference.")
        elif detected and detected[0]["Reference"] != expected_ref:
            notes.append(f"The reads match {detected[0]['Reference']} ({detected[0]['Organism']}) better "
                         f"than the configured reference {expected_ref}.")
    return notes


def screen_sample(reads: Path, out_dir: Path, sample: str, threads: int,
                  expected_ref: Optional[str] = None) -> Path:
    """STEP 5c for one sample: screen a read sample (first SCREEN_READS of R1)
    and write 02_qc/<sample>.reference_screen.tsv. Findings are logged as
    warnings; they never stop the run (the dashboard shows them)."""
    import tempfile
    out = out_dir / f"{sample}.reference_screen.tsv"
    with tempfile.TemporaryDirectory() as tmp:
        fq = Path(tmp) / "sample.fastq"
        sample_reads(reads, fq)
        rows, n = screen_reads(fq, threads)
    types = {"Organism": pl.Utf8, "Reference": pl.Utf8, "Kind": pl.Utf8, "Status": pl.Utf8,
             "Shared_reads": pl.Float64, "Evenness": pl.Float64, "Distant_evenness": pl.Float64, "Mean_identity": pl.Float64, "Reads_per_million": pl.Float64}
    pl.DataFrame(rows, schema={c: types.get(c, pl.Int64) for c in SCREEN_COLUMNS}) \
      .write_csv(out, separator="\t")
    top = next((r for r in rows if r["Status"] == "detected"), None)
    log.info(f"Reference screen ({n:,} reads): " +
             (f"best match {top['Organism']} ({top['Reference']}), spread {top['Spread']}/{SCREEN_WINDOWS}"
              if top else "no built-in genome detected"))
    for note in screen_findings(rows, expected_ref):
        log.warning(f"[{sample}] {note}")
    return out


def print_screen(rows: List[dict], n: int, label: str) -> None:
    """The screen result as a table under question 4."""
    print(f"\n      Reference screen: {n:,} reads of {label}")
    print(f"        {'Organism':<28} {'Reference':<16} {'Reads':>7} {'Shared':>7} {'Distant':>8} {'Spread':>7} {'Even':>5}  Status")
    for r in rows:
        if r["Status"] == "–" and r["Reads"] == 0 and r["Distant_reads"] == 0 and r["Shared_reads"] == 0:
            continue                          # keep the table short
        mark = "★ " if r["Status"] == "detected" else "  "
        print(f"      {mark}{r['Organism']:<28} {r['Reference']:<16} {r['Reads']:>7,} {r['Shared_reads']:>7.1f} "
              f"{r['Distant_reads']:>8,} {r['Spread']:>3}/{SCREEN_WINDOWS:<3} {r['Evenness']:>5.2f}  {r['Status']}")
    quiet = sum(1 for r in rows if r["Status"] == "–" and not (r["Reads"] or r["Distant_reads"] or r["Shared_reads"]))
    if quiet:
        print(f"        ({quiet} other references: no reads)")
    print(f"      Spread = how many of {SCREEN_WINDOWS} stretches along the genome got a close read; Even =")
    print("      spread ÷ what random placement would give (≈ 1 when the organism is really there,")
    print("      ≪ 1 when reads pile up in a few conserved genes).")


def run_reference_screen(read_source, label: str, threads: int,
                         organism_hint: str = "") -> Optional[str]:
    """Option 0 at question 4: screen, print, and return the best reference
    (None if nothing is detected)."""
    import shutil
    import tempfile
    if shutil.which("bwa") is None:
        print("      ! 'bwa' isn't on PATH. Activate the env first: source ~/miniforge3/bin/activate ngs-env")
        return None
    try:
        with tempfile.TemporaryDirectory() as tmp:
            fq = Path(tmp) / "sample.fastq"
            print(f"      Taking the first {SCREEN_READS:,} reads of {label} …", flush=True)
            if sample_reads(read_source, fq) == 0:
                print("      ! No reads could be read for the screen.")
                return None
            print(f"      Mapping them against {len(screen_panel())} genomes (first time: download + index) …", flush=True)
            rows, n = screen_reads(fq, threads)
    except ConnectionError as e:
        print(f"      ! Couldn't download the reads or references ({e}). The screen needs internet.")
        return None
    except Exception as e:
        print(f"      ! The screen failed: {e}")
        return None
    print_screen(rows, n, label)
    detected = [r for r in rows if r["Status"] == "detected"]
    hint = reference_for(organism_hint) if organism_hint else None
    for note in screen_findings(rows, hint[1] if hint else None):
        print(f"      ! {note}")
    if not detected:
        print()
        return None
    best = detected[0]
    print(f"      Best match: {best['Reference']} ({best['Organism']}), {best['Reads']:,} reads "
          f"over {best['Spread']}/{SCREEN_WINDOWS} of the genome.\n")
    return best["Reference"]


def _ask(prompt: str) -> str:
    """input() that stops cleanly when piped answers run out."""
    try:
        return input(prompt).strip()
    except EOFError:
        print("\n      ! Input ended before every question was answered.")
        sys.exit(1)

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
    print("      An accession is looked up on ENA right away (needs internet);")
    print("      a file is read once for a summary.")
    sra_info = None        # ENA's answer for an SRA accession (None: local file, or ENA unreachable)
    local_input = False
    while True:
        accession = _ask("      Accession or path: ")
        if not accession:
            print("      ! This field is required.")
            continue
        if accession.endswith(FASTQ_EXT):
            if not Path(accession).is_file():
                print(f"      ! File not found: {Path(accession).resolve()}")
                continue
            if not describe_local_fastq(Path(accession)):
                continue
            local_input = True
            break
        # WHY reject other accessions here: a genome accession (NC_…) would be
        #     downloaded as FASTA and handed to fastp as reads, which fails later.
        if not accession.upper().startswith(SRA_PREFIX):
            print(f"      ! '{accession}' isn't a read accession (SRR…, ERR…, DRR…) or a .fastq/.fq file.")
            print("        A genome accession such as NC_002549.1 goes in question 4.")
            continue
        accession = accession.upper()
        print(f"      Looking up {accession} on ENA …", flush=True)
        try:
            sra_info = lookup_sra(accession)
        except LookupError as e:
            print(f"      ! {e}. Check the spelling and try again.")
            continue
        except ConnectionError as e:
            print(f"      ! Couldn't reach ENA ({e}), so there's no preview.")
            print("        Continuing; the download in step 4 needs internet too.\n")
            break
        print_sra_info(accession, sra_info)
        break

    print(f"\n[2/4] CPU threads – this machine has {os.cpu_count()} cores")
    print("      e.g. 4   (press Enter for default 8)")
    while True:
        t = _ask("      Threads [8]: ")
        if not t:
            threads = 8
            break
        if t.isdigit() and int(t) > 0:
            threads = int(t)
            break
        print("      ! Please enter a positive whole number, e.g. 4")

    print("\n[3/4] Output folder – created if it doesn't exist")
    print("      e.g. ./NGS_run   or   ./SRR1553607_run   (press Enter for ./NGS_run)")
    base_dir = Path(_ask("      Output directory [./NGS_run]: ") or "./NGS_run")

    # Optional for SRA (looked up from the organism), required for local files
    organism_hint = (sra_info or {}).get("scientific_name") or ""
    suggested = reference_for(organism_hint) if organism_hint else None
    if suggested and suggested[1].startswith("GCF_"):
        suggested = None                 # assemblies can't be downloaded yet
    choices = reference_choices()
    print("\n[4/4] Reference genome accession (RefSeq) – optional for SRA, REQUIRED for local FASTQ")
    print("      e.g. NC_002549.1 (Zaire ebolavirus)   NC_000913.3 (E. coli K-12)")
    print("           NC_045512.2 (SARS-CoV-2)         NC_003198.1 (Salmonella enterica)")
    print("      Type  list  to see all built-in references, then pick one by number,")
    print("      or  0  to test your reads against all of them and get the best match.")
    if local_input:
        print("      Required: local files carry no organism name.")
    elif suggested:
        print(f"      Press Enter to use {suggested[1]} ({suggested[0]}), matched from the organism '{organism_hint}'.")
    elif organism_hint:
        print(f"      '{organism_hint}' has no built-in reference: type one, or  list .")
    else:
        print("      Press Enter to auto-select from the SRA organism name.")
    screen_best = None     # best match from option 0, used if Enter is pressed next
    while True:
        answer = _ask("      Reference accession: ")
        if answer.lower() in ("list", "l", "?"):
            print_reference_list()
            continue
        if answer.lower() in ("0", "all"):
            # Reads for the screen: the local R1, or the start of ENA's R1 (streamed)
            if local_input:
                source, label = Path(accession), Path(accession).name
            else:
                links = [l for l in ((sra_info or {}).get("fastq_ftp") or "").split(";") if l]
                if not links:
                    print("      ! No read files to sample (ENA wasn't reachable or lists none). Type an accession.")
                    continue
                r1_link = next((l for l in links if l.endswith("_1.fastq.gz")), links[0])
                source, label = "https://" + r1_link.replace("ftp://", ""), accession
            screen_best = run_reference_screen(source, label, threads, organism_hint) or screen_best
            if screen_best:
                print(f"      Press Enter to use {screen_best}, or type another number or accession.")
            continue
        if answer.isdigit():
            i = int(answer)
            if not 1 <= i <= len(choices):
                print(f"      ! Pick a number from 1 to {len(choices)}, or type  list .")
                continue
            answer = choices[i - 1][0]
            print(f"      → {answer} ({choices[i - 1][1]})")
        if not answer and screen_best:
            answer = screen_best
            print(f"      → {answer} (best match from the screen)")
        if not answer:
            if local_input:
                print("      ! Local FASTQ needs a reference. Type an accession, or  list  to choose.")
                continue
            if organism_hint and not suggested:
                print(f"      ! No downloadable built-in reference for '{organism_hint}'. Type one, or  list .")
                continue
            force_ref = None             # auto: download_reference() maps the organism
            break
        alias = next((k for k in REF_SETS if k.lower() == answer.lower()), None)
        answer = alias or answer.upper()
        if answer.startswith(("GCF_", "GCA_")):
            print("      ! Assembly accessions (GCF_/GCA_) can't be downloaded yet.")
            print("        Use a nucleotide accession such as NC_045512.2.")
            continue
        try:
            ref_info = check_reference(answer)
        except LookupError as e:
            print(f"      ! {e}. Check the accession, or type  list .")
            continue
        except ConnectionError as e:
            print(f"      ! Couldn't check it on NCBI ({e}); using it as typed.")
            force_ref = answer
            break
        print(f"      ✓ {answer}: {ref_info['title']} ({ref_info['length']:,} bp)")
        if suggested and suggested[1] != answer:
            print(f"      ! Your sample is '{organism_hint}', which usually uses {suggested[1]}. "
                  f"Continuing with {answer}.")
        force_ref = answer
        break
    chosen = force_ref or (suggested[1] if suggested else None)
    if chosen:
        print("      Gene annotation: " + ("yes (SnpEff " + SNPEFF_DB[chosen] + ")" if chosen in SNPEFF_DB
              else "skipped – no SnpEff database for this reference (variants are still reported)"))
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
        mate_re = MATE_RE
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

    # ----- STEP 5b–5d · Sample checks -----
    # WHY warn and continue: these checks inform you; the variant calling can
    #     still run. (Snakemake/Nextflow stop on a failing step instead.)
    if HOST_INDEX:
        try:
            r1_trim, r2_trim = remove_host(r1_trim, r2_trim, HOST_INDEX, dirs["qc"], sample, threads)
        except Exception as e:
            log.warning(f"Host removal failed ({e}); continuing with all reads.")
    else:
        log.info("Host removal: off (set NGS_HOST_INDEX to a bwa-indexed host genome to turn it on)")
    if os.environ.get("NGS_SCREEN", "1") not in ("0", "false", "False"):
        try:
            screen_sample(r1_trim, dirs["qc"], sample, threads, expected_ref=ref_fasta.stem)
        except Exception as e:
            log.warning(f"Reference screen failed ({e}); continuing.")
    if KRAKEN_DB:
        try:
            classify_reads(r1_trim, r2_trim, KRAKEN_DB, dirs["qc"], sample, threads)
        except Exception as e:
            log.warning(f"Kraken2 classification failed ({e}); continuing.")
    else:
        log.info("Kraken2 classification: off (set NGS_KRAKEN_DB to a Kraken2 database to turn it on)")

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