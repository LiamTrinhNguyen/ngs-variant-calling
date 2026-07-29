#!/usr/bin/env python3
"""
NGS_Pipeline_v2.6 – Full automated variant calling pipeline
Author : Liam TrinhNguyen
Email  : LiamTrinhNguyen@gmail.com
"""

__author__  = "Liam TrinhNguyen"
__email__   = "LiamTrinhNguyen@gmail.com"
__version__ = "NGS_Pipeline_v2.6"

import os
import sys
import logging
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Optional, List, Tuple

import polars as pl

# ---------------------------------------------------------------------------
# Optional / auto-install friendly imports
# ---------------------------------------------------------------------------
try:
    from tqdm import tqdm
except ImportError:
    print("Installing tqdm...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "tqdm", "--quiet"])
    from tqdm import tqdm

try:
    import requests
except ImportError:
    requests = None

try:
    from Bio import Entrez, SeqIO
    Entrez.email = "LiamTrinhNguyen@gmail.com"
except ImportError:
    Entrez = None
    SeqIO = None

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s | %(levelname)-8s | %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
log = logging.getLogger("NGS_Pipeline")

# ---------------------------------------------------------------------------
# 1. SETUP – organised directory structure
# ---------------------------------------------------------------------------
def setup_directories(base_dir: Path) -> dict:
    """Create a reproducible project layout."""
    dirs = {
        "base":      base_dir,
        "raw":       base_dir / "00_raw",
        "ref":       base_dir / "01_reference",
        "qc":        base_dir / "02_qc",
        "aligned":   base_dir / "03_aligned",
        "variants":  base_dir / "04_variants",
        "annotated": base_dir / "05_annotated",
        "reports":   base_dir / "06_reports",
        "logs":      base_dir / "logs",
    }
    for d in dirs.values():
        d.mkdir(parents=True, exist_ok=True)
    log.info(f"Project directories created under {base_dir}")
    return dirs

# ---------------------------------------------------------------------------
# Helper: run external tools safely
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
# Accession helpers
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

    return organism, fastq_paths

# ---------------------------------------------------------------------------
# 2. REFERENCE – organism-driven reference genome
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
    if not fasta.exists():
        download_ncbi_fasta(acc, fasta)
    else:
        log.info(f"Reference already present: {fasta}")

    return fasta

# ---------------------------------------------------------------------------
# 3. INDEXING – BWA + samtools
# ---------------------------------------------------------------------------
def index_reference(fasta: Path, log_dir: Path) -> None:
    log.info(f"Building BWA index for {fasta.name}")
    run_cmd(["bwa", "index", str(fasta)], log_dir / "bwa_index.log")

    log.info("Building samtools faidx")
    run_cmd(["samtools", "faidx", str(fasta)], log_dir / "samtools_faidx.log")

    dict_file = fasta.with_suffix(".dict")
    if not dict_file.exists():
        run_cmd(
            ["samtools", "dict", str(fasta), "-o", str(dict_file)],
            log_dir / "samtools_dict.log",
        )

# ---------------------------------------------------------------------------
# 5. fastp – QC + trimming + base correction
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
        "--correction",
        "--thread", str(threads),
        "--qualified_quality_phred", "20",
        "--length_required", "30",
    ]
    if r2:
        cmd.extend(["-I", str(r2), "-O", str(out_r2), "--detect_adapter_for_pe"])
    else:
        cmd.append("--detect_adapter_for_se")

    run_cmd(cmd, out_dir / f"{sample}.fastp.log")
    return out_r1, out_r2

# ---------------------------------------------------------------------------
# 6. BWA-MEM + Samtools – alignment → sorted BAM
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
        "-R", f"@RG\\tID:{sample}\\tSM:{sample}\\tPL:ILLUMINA",
        str(ref), str(r1)
    ]
    if r2:
        bwa_cmd.append(str(r2))

    log.info("Running BWA-MEM")
    with open(sam, "w") as fh:
        subprocess.run(bwa_cmd, stdout=fh, check=True)

    log.info("Sorting & indexing BAM")
    run_cmd(["samtools", "sort", "-@", str(threads), "-o", str(bam), str(sam)])
    run_cmd(["samtools", "index", str(bam)])
    sam.unlink(missing_ok=True)
    return bam

# ---------------------------------------------------------------------------
# 7. BCFtools – variant calling
# ---------------------------------------------------------------------------
def call_variants(ref: Path, bam: Path, out_dir: Path, sample: str) -> Path:
    vcf = out_dir / f"{sample}.raw.vcf.gz"
    log.info("Calling variants with bcftools mpileup + call")

    cmd = (
        f"bcftools mpileup -f {ref} -Ou {bam} | "
        f"bcftools call -mv -Oz -o {vcf}"
    )
    run_cmd(["bash", "-c", cmd], out_dir / f"{sample}.bcftools.log")
    run_cmd(["bcftools", "index", str(vcf)])
    return vcf

# ---------------------------------------------------------------------------
# 8. SnpEff – functional annotation
# ---------------------------------------------------------------------------
def annotate_variants(vcf: Path, genome: str, out_dir: Path, sample: str) -> Path:
    """
    genome should be a SnpEff database name.
    Example for Ebola: you may need to build a custom database or use a close relative.
    """
    ann_vcf = out_dir / f"{sample}.ann.vcf.gz"
    log.info(f"Annotating with SnpEff ({genome})")

    # Write uncompressed first
    uncompressed = ann_vcf.with_suffix("")
    cmd = ["snpEff", "-v", genome, str(vcf)]
    with open(uncompressed, "w") as fh:
        subprocess.run(cmd, stdout=fh, check=True)

    run_cmd(["bgzip", "-f", str(uncompressed)])
    run_cmd(["tabix", "-p", "vcf", str(ann_vcf)])
    return ann_vcf

# ---------------------------------------------------------------------------
# 9. Polars – high-confidence filtering & reporting
# ---------------------------------------------------------------------------
def filter_and_report(vcf: Path, reports_dir: Path, sample: str, min_qual: float = 30.0) -> Path:
    """
    Parse VCF (annotated or raw) with Polars and produce a clean TSV report.
    Handles both SnpEff-annotated VCFs (with INFO/ANN) and raw bcftools VCFs.
    """
    log.info(f"Parsing VCF with Polars: {vcf.name}")

    tmp_tsv = reports_dir / f"{sample}.tmp.tsv"

    # First try the annotated format
    try:
        run_cmd([
            "bcftools", "query",
            "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/ANN\n",
            str(vcf),
            "-o", str(tmp_tsv),
        ], check=True)
        has_ann = True
    except RuntimeError:
        # Fall back to raw VCF format (no ANN field)
        log.warning("No INFO/ANN field found – treating as raw (unannotated) VCF")
        run_cmd([
            "bcftools", "query",
            "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\n",
            str(vcf),
            "-o", str(tmp_tsv),
        ])
        has_ann = False

    if has_ann:
        df = pl.read_csv(
            tmp_tsv,
            separator="\t",
            has_header=False,
            new_columns=["CHROM", "POS", "REF", "ALT", "QUAL", "ANN"],
            ignore_errors=True,
        )
        df = df.filter(pl.col("QUAL").cast(pl.Float64, strict=False) >= min_qual)

        df = df.with_columns([
            pl.col("ANN").str.split("|").list.get(1).alias("Impact"),
            pl.col("ANN").str.split("|").list.get(3).alias("Gene"),
            pl.col("ANN").str.split("|").list.get(9).alias("HGVS_c"),
            pl.col("ANN").str.split("|").list.get(10).alias("HGVS_p"),
        ]).select([
            "CHROM", "POS", "REF", "ALT", "QUAL",
            "Impact", "Gene", "HGVS_c", "HGVS_p"
        ])
    else:
        df = pl.read_csv(
            tmp_tsv,
            separator="\t",
            has_header=False,
            new_columns=["CHROM", "POS", "REF", "ALT", "QUAL"],
            ignore_errors=True,
        )
        df = df.filter(pl.col("QUAL").cast(pl.Float64, strict=False) >= min_qual)
        # Add empty annotation columns for consistent output
        df = df.with_columns([
            pl.lit(None).alias("Impact"),
            pl.lit(None).alias("Gene"),
            pl.lit(None).alias("HGVS_c"),
            pl.lit(None).alias("HGVS_p"),
        ])

    report = reports_dir / f"{sample}.high_confidence_variants.tsv"
    df.write_csv(report, separator="\t")
    log.info(f"High-confidence report written → {report}  ({df.height} variants)")
    tmp_tsv.unlink(missing_ok=True)
    return report

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    print(f"\n{'='*65}")
    print(f"  NGS Pipeline {__version__}")
    print(f"  Author: {__author__}  <{__email__}>")
    print(f"{'='*65}\n")

    accession = input("Accession (SRR… / NM_… / local R1.fastq.gz path): ").strip()
    threads   = int(input("Threads to use [8]: ") or 8)
    base_dir  = Path(input("Project output directory [./NGS_run]: ").strip() or "./NGS_run")

    # Optional: force a specific reference
    force_ref = input("Force reference accession? (leave empty for auto) : ").strip() or None

    dirs = setup_directories(base_dir)

    # ----- 4. Data Input -----
    if accession.endswith((".fastq", ".fastq.gz", ".fq", ".fq.gz")):
        r1 = Path(accession)
        if not r1.exists():
            raise FileNotFoundError(f"Local FASTQ not found: {r1}")
        # simple paired-end detection
        r2_candidates = [
            r1.parent / r1.name.replace("_R1", "_R2").replace("_1", "_2"),
            r1.parent / r1.name.replace("R1", "R2"),
        ]
        r2 = next((p for p in r2_candidates if p.exists()), None)
        organism = "unknown"
        sample = r1.stem.split(".")[0]
        log.info(f"Using local FASTQ: {r1}" + (f" + {r2}" if r2 else " (single-end)"))
    else:
        organism, fastqs = process_accession(accession, dirs)
        sample = accession
        r1 = fastqs[0]
        r2 = fastqs[1] if len(fastqs) > 1 else None

    # ----- 2. Reference -----
    ref_fasta = download_reference(organism, dirs["ref"], force_accession=force_ref)

    # ----- 3. Indexing -----
    index_reference(ref_fasta, dirs["logs"])

    # ----- 5. fastp -----
    r1_trim, r2_trim = run_fastp(r1, r2, dirs["qc"], sample, threads)

    # ----- 6. Alignment -----
    bam = align_and_sort(ref_fasta, r1_trim, r2_trim, dirs["aligned"], sample, threads)

    # ----- 7. Variant calling -----
    raw_vcf = call_variants(ref_fasta, bam, dirs["variants"], sample)

    # ----- 8. Annotation -----
    # Note: For Zaire ebolavirus you may need to build a custom SnpEff database
    # or temporarily skip annotation.
    snpeff_genome = "Escherichia_coli_K_12_substr_MG1655"   # change if you have a better DB
    try:
        ann_vcf = annotate_variants(raw_vcf, snpeff_genome, dirs["annotated"], sample)
    except Exception as e:
        log.warning(f"SnpEff annotation skipped ({e}). Using raw VCF for reporting.")
        ann_vcf = raw_vcf

    # ----- 9. Filtering & reporting -----
    report = filter_and_report(ann_vcf, dirs["reports"], sample)

    print("\n" + "=" * 65)
    print("PIPELINE COMPLETED SUCCESSFULLY")
    print(f"  BAM           : {bam}")
    print(f"  Annotated VCF : {ann_vcf}")
    print(f"  Final report  : {report}")
    print("=" * 65)

if __name__ == "__main__":
    main()