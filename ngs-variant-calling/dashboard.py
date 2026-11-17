#!/usr/bin/env python3
"""
STEP 10 · DASHBOARD – turn a run's output folder into easy-to-read results.

WHAT: reads every sample in a run folder (00_raw … 06_reports) and writes
      07_dashboard/
        dashboard.html          one interactive page: QC checks, coverage,
                                variants along the genome, filters, variant table,
                                and both CSVs as searchable "Data tables"
        sample_summary.csv      one row per sample   (for Power BI / Excel)
        all_samples_variants.csv  every kept (PASS) variant, with a Sample column
        all_samples_calls.csv   every call with its FILTER reason (audit trail)
        r_statistics.csv        the R statistics (step 9b), one row per test
        variant_windows.csv     variants per genome window + hotspot test (from R)
WHY:  the TSV/VCF/BAM files are hard to read. The page shows at a glance
      whether a sample is good and what was found. The CSVs can be read in the
      page (no download needed) and also load straight into Power BI, Excel or R.

The HTML is self-contained (data + code in one file, no internet needed), so it
can be emailed or opened from a shared drive.

Usage:
  python3 dashboard.py <run folder> [--min-qual 30] [--min-dp 15] [--min-af 0.8]
View (Codespaces):
  python3 -m http.server 8000 --directory <run folder>/07_dashboard &
  "$BROWSER" http://localhost:8000/dashboard.html
  (also runs automatically at the end of the interactive script, Snakemake and Nextflow)
"""

import argparse
import csv
import json
import math
import re
import subprocess
from datetime import datetime
from pathlib import Path


# ---------------------------------------------------------------------------
# Helpers – read the files each step produced
# ---------------------------------------------------------------------------
def _run(cmd):
    """Run a tool and return its stdout (empty string if it fails)."""
    r = subprocess.run(cmd, capture_output=True, text=True)
    return r.stdout if r.returncode == 0 else ""


def _reference_info(run_dir: Path):
    """Reference name, description and contig lengths (from the .fai index)."""
    fai = next((run_dir / "01_reference").glob("*.fai"), None)
    if fai is None:
        return {"accession": "unknown", "description": "", "contigs": []}
    fasta = fai.with_suffix("")
    desc = ""
    if fasta.exists():
        with open(fasta) as fh:
            desc = fh.readline().lstrip(">").strip()
    contigs = []
    for line in fai.read_text().splitlines():
        name, length = line.split("\t")[:2]
        contigs.append({"name": name, "length": int(length)})
    return {"accession": fasta.stem, "description": desc, "contigs": contigs}


def _fastp(json_path: Path):
    """Read counts and quality before/after trimming."""
    if not json_path.exists():
        return None
    d = json.loads(json_path.read_text())
    b, a = d["summary"]["before_filtering"], d["summary"]["after_filtering"]
    f = d.get("filtering_result", {})
    return {
        "reads_raw": b["total_reads"],
        "reads_pass": a["total_reads"],
        "q30_before": b["q30_rate"],
        "q30_after": a["q30_rate"],
        "gc": a["gc_content"],
        "low_quality": f.get("low_quality_reads", 0),
        "too_short": f.get("too_short_reads", 0),
    }


def _settings(run_dir: Path, name: str):
    """How the variants were called: read from 04_variants/<sample>.calling_settings.json,
    which the calling step writes with the exact options it used.
    WHY that file: bcftools mpileup doesn't record its options in the VCF header,
         so the file is the only reliable record. Runs made before it existed
         show "unknown" rather than a guess."""
    f = run_dir / "04_variants" / f"{name}.calling_settings.json"
    if f.exists():
        return json.loads(f.read_text())
    # Older run: only haploid / normalized can be read from the VCF header
    vcf = run_dir / "04_variants" / f"{name}.raw.vcf.gz"
    head = _run(["bcftools", "view", "-h", str(vcf)]) if vcf.exists() else ""
    return {
        "min_mapq": None, "min_baseq": None, "max_depth": None,
        "haploid": "--ploidy 1" in head,
        "normalized": "##bcftools_normCommand" in head,
        "allele_depths": None,
    }


def _alignment(bam: Path, s):
    """Mapping rate (samtools flagstat) and coverage (samtools coverage).
    Coverage uses the SAME read/base quality cut-offs as variant calling.
    WHY: the depth shown then matches the depth the filter uses."""
    if not bam.exists():
        return None, None
    out = {}
    fs = _run(["samtools", "flagstat", "-O", "json", str(bam)])
    if fs:
        q = json.loads(fs)["QC-passed reads"]
        out = {
            "reads": q["primary"],
            "mapped": q["primary mapped"],
            "mapped_pct": q["primary mapped %"],
            "proper_pct": q.get("properly paired %"),
        }
    cov = {"meandepth": 0.0, "covered_pct": 0.0, "meanmapq": 0.0}
    rows = [l.split("\t") for l in _run(["samtools", "coverage", "-q", str(s["min_mapq"] or 0),
                                         "-Q", str(s["min_baseq"] or 0), str(bam)]).splitlines()[1:]]
    if rows:
        # Combine contigs, weighted by length (bacteria can have plasmids)
        total = sum(int(r[2]) - int(r[1]) + 1 for r in rows)
        covered = sum(int(r[4]) for r in rows)
        depth_sum = sum(float(r[6]) * (int(r[2]) - int(r[1]) + 1) for r in rows)
        mapq_sum = sum(float(r[8]) * (int(r[2]) - int(r[1]) + 1) for r in rows)
        cov = {
            "meandepth": depth_sum / total,
            "covered_pct": 100 * covered / total,
            "meanmapq": mapq_sum / total,
        }
    return out, cov


def _depth_profile(bam: Path, ref, s, n_bins: int = 600):
    """Mean good-quality depth in ~600 equal bins along the genome (contigs joined)."""
    genome_len = sum(c["length"] for c in ref["contigs"])
    if not bam.exists() or genome_len == 0:
        return None
    n_bins = min(n_bins, genome_len)
    bin_size = math.ceil(genome_len / n_bins)
    offsets, off = {}, 0
    for c in ref["contigs"]:
        offsets[c["name"]] = off
        off += c["length"]
    sums = [0] * n_bins
    counts = [0] * n_bins
    # -q/-Q: same cut-offs as variant calling, so the plot matches the depth filter
    proc = subprocess.Popen(["samtools", "depth", "-a", "-q", str(s["min_baseq"] or 0),
                             "-Q", str(s["min_mapq"] or 0), str(bam)],
                            stdout=subprocess.PIPE, text=True)
    for line in proc.stdout:
        chrom, pos, depth = line.split("\t")
        i = min((offsets.get(chrom, 0) + int(pos) - 1) // bin_size, n_bins - 1)
        sums[i] += int(depth)
        counts[i] += 1
    proc.wait()
    return {
        "bin_size": bin_size,
        "depth": [round(s / c, 2) if c else 0 for s, c in zip(sums, counts)],
    }


def _num(v, cast=float, default=None):
    try:
        return cast(v)
    except (TypeError, ValueError):
        return default


def _calls(run_dir: Path, name: str, offsets, t):
    """Every call with its FILTER reasons, from 06_reports/<sample>.all_calls.tsv.
    WHY that file: it holds the pipeline's own decisions, so the charts show
         exactly what the filter did instead of re-computing it here."""
    tsv = run_dir / "06_reports" / f"{name}.all_calls.tsv"
    calls = []
    if tsv.exists():
        for r in _report_rows(tsv):
            f = r.get("FILTER") or ""
            calls.append({
                "chrom": r["CHROM"], "pos": int(r["POS"]), "gpos": offsets.get(r["CHROM"], 0) + int(r["POS"]),
                "ref": r["REF"], "alt": r["ALT"], "qual": round(_num(r["QUAL"], default=0.0), 1),
                "dp": _num(r["DP"], int, 0), "af": _num(r.get("AF")),
                "filter": f, "kept": f == "PASS",
            })
        return calls
    # Older runs (before flagging was added): re-compute from the raw VCF
    vcf = run_dir / "04_variants" / f"{name}.raw.vcf.gz"
    if not vcf.exists():
        return []
    for line in _run(["bcftools", "query", "-f", "%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/DP\n", str(vcf)]).splitlines():
        chrom, pos, ref, alt, qual, dp = line.split("\t")
        q, d = _num(qual, default=0.0), _num(dp, int, 0)
        reasons = [n for n, bad in (("LowQual", q < t["min_qual"]), ("LowDepth", d < t["min_dp"])) if bad]
        calls.append({
            "chrom": chrom, "pos": int(pos), "gpos": offsets.get(chrom, 0) + int(pos),
            "ref": ref, "alt": alt, "qual": round(q, 1), "dp": d, "af": None,
            "filter": ";".join(reasons) or "PASS", "kept": not reasons,
        })
    return calls


def _r_stats(run_dir: Path, name: str):
    """Statistics from R (step 9b): 06_reports/<sample>.r_stats.tsv and .r_windows.tsv.
    WHY Polars here: R writes plain TSV; Polars reads it with the right types
         (numbers stay numbers, the TRUE/FALSE Hotspot column becomes a boolean)
         and turns it into rows for the page. Returns None for runs made
         before the R step existed, and the page then hides the section."""
    import polars as pl   # WHY here: only this optional section needs it
    stats_f = run_dir / "06_reports" / f"{name}.r_stats.tsv"
    win_f = run_dir / "06_reports" / f"{name}.r_windows.tsv"
    if not stats_f.exists():
        return None
    stats = pl.read_csv(stats_f, separator="\t", infer_schema_length=0)   # all text: values are labels like "5.57 (3.55–9.09)"
    windows = pl.read_csv(win_f, separator="\t") if win_f.exists() else pl.DataFrame()
    if "Hotspot" in windows.columns and windows["Hotspot"].dtype == pl.Utf8:
        windows = windows.with_columns(pl.col("Hotspot") == "TRUE")
    return {"stats": stats.to_dicts(), "windows": windows.to_dicts()}


def _report_rows(tsv: Path):
    with open(tsv, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------------------
# Collect all samples in a run folder
# ---------------------------------------------------------------------------
def collect(run_dir: Path, min_qual: float, min_dp: int, min_af: float, only=None):
    ref = _reference_info(run_dir)
    offsets, off = {}, 0
    for c in ref["contigs"]:
        offsets[c["name"]] = off
        off += c["length"]
    t = {"min_qual": min_qual, "min_dp": min_dp, "min_af": min_af}

    samples = []
    # WHY start from the reports: a sample is only "done" once it has a report
    for tsv in sorted((run_dir / "06_reports").glob("*.high_confidence_variants.tsv")):
        name = tsv.name.replace(".high_confidence_variants.tsv", "")
        if only is not None and name not in only:
            continue   # e.g. a sample removed from the Snakemake config
        bam = run_dir / "03_aligned" / f"{name}.sorted.bam"
        settings = _settings(run_dir, name)
        align, cov = _alignment(bam, settings)
        variants = _report_rows(tsv)
        samples.append({
            "name": name,
            "settings": settings,
            "qc": _fastp(run_dir / "02_qc" / f"{name}.fastp.json"),
            "align": align,
            "cov": cov,
            "profile": _depth_profile(bam, ref, settings),
            "calls": _calls(run_dir, name, offsets, t),
            "variants": variants,
            "annotated": any(v.get("Effect") for v in variants),
            "r": _r_stats(run_dir, name),
        })
    return {
        "run": str(run_dir.resolve()),
        "generated": datetime.now().strftime("%Y-%m-%d %H:%M"),
        "reference": ref,
        "offsets": offsets,
        "thresholds": t,
        "samples": samples,
    }


# ---------------------------------------------------------------------------
# CSV tables – shown in the page's "Data tables" section AND saved as files
#   for Power BI / Excel / R
#   WHY long format (one row per variant, with a Sample column): Power BI and
#   Excel pivot tables can then slice by sample, gene or impact directly.
# ---------------------------------------------------------------------------
def build_tables(data):
    """Rows for the three CSV files.
    WHY one function: the same rows are written to disk AND shown as tables in
         the dashboard, so the page and the files can never disagree."""
    t = data["thresholds"]
    summary_cols = ["Sample", "Reference", "Reads_raw", "Reads_after_QC", "Pct_reads_kept_QC",
                    "Q30_after_QC", "Mapped_reads", "Pct_mapped", "Mean_depth",
                    "Pct_genome_covered", "Raw_calls", "Flagged_LowQual", "Flagged_LowDepth",
                    "Flagged_LowAF", "Mixed_calls_AF_0.2_to_min", "High_confidence_variants", "Annotated"]
    summary_rows = []
    for s in data["samples"]:
        qc, al, cv = s["qc"] or {}, s["align"] or {}, s["cov"] or {}
        calls = s["calls"]
        # WHY count reasons: shows which filter removed most calls for this sample
        reason = lambda r: sum(r in c["filter"] for c in calls)
        # WHY "mixed calls": several calls with 20–80% of reads suggest a mixed
        #     sample or contamination, even if the filter removed them
        mixed = sum(c["af"] is not None and 0.2 <= c["af"] < t["min_af"] for c in calls)
        summary_rows.append([
            s["name"], data["reference"]["accession"],
            qc.get("reads_raw"), qc.get("reads_pass"),
            round(100 * qc["reads_pass"] / qc["reads_raw"], 2) if qc.get("reads_raw") else None,
            round(100 * qc["q30_after"], 2) if qc else None,
            al.get("mapped"), al.get("mapped_pct"),
            round(cv.get("meandepth", 0), 2), round(cv.get("covered_pct", 0), 2),
            len(calls), reason("LowQual"), reason("LowDepth"), reason("LowAF"), mixed,
            len(s["variants"]), s["annotated"],
        ])
    var_cols = ["Sample", "CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AF",
                "Effect", "Impact", "Gene", "HGVS_c", "HGVS_p"]
    var_rows = [[s["name"]] + [v.get(c, "") for c in var_cols[1:]]
                for s in data["samples"] for v in s["variants"]]
    # WHY a table of ALL calls: the audit trail – what was excluded and why
    call_cols = ["Sample", "CHROM", "POS", "REF", "ALT", "QUAL", "DP", "AF", "FILTER"]
    call_rows = [[s["name"], c["chrom"], c["pos"], c["ref"], c["alt"], c["qual"], c["dp"], c["af"], c["filter"]]
                 for s in data["samples"] for c in s["calls"]]
    # WHY two tables for the R step: the statistics are one row per test, the
    #     windows one row per region; both are long format with a Sample column
    r_cols = ["Sample", "Metric", "Value", "Interpretation"]
    r_rows = [[s["name"], r["Metric"], r["Value"], r["Interpretation"]]
              for s in data["samples"] if s.get("r") for r in s["r"]["stats"]]
    w_cols = ["Sample", "Window_start", "Window_end", "Variants", "Expected", "P_value", "Q_value", "Hotspot"]
    w_rows = [[s["name"]] + [w.get(c) for c in w_cols[1:]]
              for s in data["samples"] if s.get("r") for w in s["r"]["windows"]]
    return {
        "summary":  {"file": "sample_summary.csv",       "label": "Sample summary",           "columns": summary_cols, "rows": summary_rows},
        "variants": {"file": "all_samples_variants.csv", "label": "Trusted variants (PASS)",  "columns": var_cols,     "rows": var_rows},
        "calls":    {"file": "all_samples_calls.csv",    "label": "All calls + filter reason", "columns": call_cols,    "rows": call_rows},
        "rstats":   {"file": "r_statistics.csv",         "label": "R statistics",              "columns": r_cols,       "rows": r_rows},
        "windows":  {"file": "variant_windows.csv",      "label": "Hotspot windows (R)",       "columns": w_cols,       "rows": w_rows},
    }


def write_csvs(tables, out_dir: Path):
    paths = []
    for t in tables.values():
        p = out_dir / t["file"]
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(t["columns"])
            w.writerows(t["rows"])
        paths.append(p)
    return paths


# ---------------------------------------------------------------------------
# Build everything
# ---------------------------------------------------------------------------
def build_dashboard(run_dir, min_qual: float = 30.0, min_dp: int = 15, min_af: float = 0.8,
                    samples=None) -> Path:
    """Write 07_dashboard/ (HTML + CSVs) for a run folder; returns the HTML path.

    samples: optional list of sample names to include (default: every sample
             with a report in 06_reports/).
    """
    run_dir = Path(run_dir)
    out_dir = run_dir / "07_dashboard"
    out_dir.mkdir(parents=True, exist_ok=True)
    data = collect(run_dir, min_qual, min_dp, min_af, only=set(samples) if samples else None)
    data["tables"] = build_tables(data)       # shown in the page's "Data tables" section
    write_csvs(data["tables"], out_dir)        # the same rows, saved for Power BI / Excel
    # WHY escape "</": stops sample text from closing the <script> tag early
    payload = json.dumps(data, separators=(",", ":")).replace("</", "<\\/")
    html = out_dir / "dashboard.html"
    html.write_text(HTML_TEMPLATE.replace("__DATA__", payload), encoding="utf-8")
    return html


# ---------------------------------------------------------------------------
# The page (HTML + CSS + JS in one string; data is inserted at __DATA__)
# ---------------------------------------------------------------------------
HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Variant calling results</title>
<style>
  :root {
    color-scheme: light;
    --bg: #f4f4f1; --surface: #fcfcfb; --border: #e3e2dd; --grid: #e8e7e3;
    --text: #0b0b0b; --text-2: #52514e; --muted: #7c7b75;
    --accent: #2a78d6; --accent-wash: rgba(42,120,214,.10);
    --good: #0ca30c; --warning: #fab219;
  }
  @media (prefers-color-scheme: dark) {
    :root {
      color-scheme: dark;
      --bg: #121211; --surface: #1a1a19; --border: #2f2f2d; --grid: #2b2b29;
      --text: #f0efec; --text-2: #c3c2b7; --muted: #8f8e86;
      --accent: #3987e5; --accent-wash: rgba(57,135,229,.14);
    }
  }
  * { box-sizing: border-box; }
  body { margin: 0; background: var(--bg); color: var(--text);
         font: 14px/1.5 -apple-system, "Segoe UI", Helvetica, Arial, sans-serif; }
  .wrap { max-width: 1200px; margin: 0 auto; padding: 24px 16px 48px; }
  header h1 { margin: 0 0 4px; font-size: 22px; }
  header p { margin: 0; color: var(--text-2); }
  .filters { display: flex; flex-wrap: wrap; gap: 12px; align-items: center; margin: 20px 0 16px; }
  .filters label { color: var(--text-2); font-size: 13px; }
  select, input[type=search] { font: inherit; padding: 6px 10px; border-radius: 8px;
         border: 1px solid var(--border); background: var(--surface); color: var(--text); }
  button { font: inherit; padding: 6px 12px; border-radius: 8px; border: 1px solid var(--border);
           background: var(--surface); color: var(--text); cursor: pointer; }
  button:hover { border-color: var(--accent); }
  .card { background: var(--surface); border: 1px solid var(--border); border-radius: 12px;
          padding: 16px 18px; margin-bottom: 16px; }
  .card h2 { margin: 0 0 2px; font-size: 15px; }
  .card .sub { margin: 0 0 12px; color: var(--text-2); font-size: 13px; }
  .tiles { display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 12px; margin-bottom: 16px; }
  .tile { background: var(--surface); border: 1px solid var(--border); border-radius: 12px; padding: 12px 14px; }
  .tile .label { color: var(--text-2); font-size: 12px; }
  .tile .value { font-size: 24px; font-weight: 600; margin-top: 2px; }
  .tile .note { color: var(--muted); font-size: 12px; }
  .grid2 { display: grid; grid-template-columns: 1fr 1fr; gap: 16px; }
  @media (max-width: 820px) { .grid2 { grid-template-columns: 1fr; } }
  .grid2 .card { margin-bottom: 0; }
  #rTable td { white-space: normal; vertical-align: top; }
  #rTable td.why { color: var(--text-2); font-size: 12px; }
  .checks { list-style: none; margin: 0; padding: 0; display: grid; gap: 8px; }
  .checks li { display: flex; gap: 10px; align-items: flex-start; }
  .badge { flex: none; display: inline-flex; align-items: center; gap: 4px; font-size: 12px; font-weight: 600;
           padding: 1px 8px; border-radius: 999px; border: 1px solid var(--border); white-space: nowrap; }
  .badge .dot { width: 8px; height: 8px; border-radius: 50%; }
  .checks .why { color: var(--text-2); font-size: 13px; }
  .table-wrap { overflow-x: auto; }
  table { width: 100%; border-collapse: collapse; font-size: 13px; }
  th, td { text-align: left; padding: 6px 10px; border-bottom: 1px solid var(--grid); white-space: nowrap; }
  th { color: var(--text-2); font-weight: 600; cursor: pointer; user-select: none; position: sticky; top: 0; background: var(--surface); }
  th.num, td.num { text-align: right; font-variant-numeric: tabular-nums; }
  tr.clickable { cursor: pointer; }
  tr.selected td { background: var(--accent-wash); }
  .scroll { max-height: 460px; overflow: auto; }
  .legend { display: flex; flex-wrap: wrap; gap: 14px; font-size: 12px; color: var(--text-2); margin-bottom: 6px; }
  .legend span { display: inline-flex; align-items: center; gap: 6px; }
  .legend i { width: 10px; height: 10px; border-radius: 50%; display: inline-block; }
  .chart { width: 100%; }
  .chart svg { display: block; }
  .empty { color: var(--text-2); font-size: 13px; padding: 12px 0; }
  #tip { position: fixed; pointer-events: none; z-index: 10; background: var(--surface); color: var(--text);
         border: 1px solid var(--border); border-radius: 8px; padding: 6px 10px; font-size: 12px;
         box-shadow: 0 4px 16px rgba(0,0,0,.12); display: none; max-width: 260px; }
  #tip b { font-weight: 600; }
  .btnlink { color: var(--accent); text-decoration: none; font-weight: 600; font-size: 13px; }
  .btnlink:hover { text-decoration: underline; }
  .tabs { display: inline-flex; border: 1px solid var(--border); border-radius: 8px; overflow: hidden; }
  .tabs button { border: 0; border-radius: 0; background: var(--surface); color: var(--text-2); }
  .tabs button + button { border-left: 1px solid var(--border); }
  .tabs button[aria-selected="true"] { background: var(--accent-wash); color: var(--text); font-weight: 600; }
  footer { color: var(--muted); font-size: 12px; margin-top: 24px; }
  code { font-size: 12px; }
</style>
</head>
<body>
<div class="wrap">
  <header>
    <h1>Variant calling results</h1>
    <p id="meta"></p>
  </header>

  <div class="filters">
    <label for="sampleSel">Sample</label>
    <select id="sampleSel"></select>
    <span style="flex:1"></span>
    <a class="btnlink" href="#dataTables">View data tables (CSV) ↓</a>
  </div>

  <div class="card" id="samplesCard">
    <h2>All samples</h2>
    <p class="sub">Click a row to view that sample below.</p>
    <div class="table-wrap"><table id="samplesTable"></table></div>
  </div>

  <div class="tiles" id="tiles"></div>

  <div class="grid2" style="margin-bottom:16px">
    <div class="card">
      <h2>Quality checks</h2>
      <p class="sub">Is this sample good enough to trust the variants?</p>
      <ul class="checks" id="checks"></ul>
    </div>
    <div class="card">
      <h2>From raw calls to trusted variants</h2>
      <p class="sub" id="funnelSub"></p>
      <div class="chart" id="funnel"></div>
    </div>
  </div>

  <div class="card">
    <h2>Read depth along the genome</h2>
    <p class="sub" id="covSub"></p>
    <div class="chart" id="coverage"></div>
  </div>

  <div class="grid2" style="margin-bottom:16px">
    <div class="card">
      <h2>Which calls pass the filters</h2>
      <p class="sub" id="scatterSub"></p>
      <div class="legend"><span><i id="lgKept"></i>Kept</span><span><i id="lgDrop"></i>Filtered out</span></div>
      <div class="chart" id="scatter"></div>
    </div>
    <div class="card">
      <h2 id="typesTitle">Types of change</h2>
      <p class="sub" id="typesSub"></p>
      <div class="chart" id="types"></div>
    </div>
  </div>

  <div class="grid2" style="margin-bottom:16px">
    <div class="card">
      <h2>Allele frequency of every call</h2>
      <p class="sub" id="afSub"></p>
      <div class="chart" id="afhist"></div>
    </div>
    <div class="card">
      <h2>Why calls were filtered out</h2>
      <p class="sub" id="reasonSub"></p>
      <div class="chart" id="reasons"></div>
    </div>
  </div>

  <div class="card" id="rCard">
    <h2>Statistics from R</h2>
    <p class="sub">Computed by <code>R/variant_stats.R</code> (step 9b) from this sample's calls. Polars builds the tables; R runs the tests.</p>
    <div class="grid2" style="margin-top:8px">
      <div class="table-wrap"><table id="rTable"></table></div>
      <div>
        <p class="sub" id="winSub"></p>
        <div class="legend"><span><i id="lgHot"></i>Hotspot (BH q &lt; 0.05)</span><span><i id="lgWin"></i>Other windows</span><span><i id="lgExp"></i>Expected if even</span></div>
        <div class="chart" id="windows"></div>
      </div>
    </div>
  </div>

  <div class="card">
    <h2>How the variants were called</h2>
    <p class="sub">Read from the VCF header, so this is exactly what ran for the selected sample.</p>
    <ul class="checks" id="settings"></ul>
  </div>

  <div class="card">
    <h2>High-confidence variants</h2>
    <p class="sub" id="varSub"></p>
    <div class="filters" style="margin:0 0 10px">
      <input type="search" id="varSearch" placeholder="Search position, gene, effect…" aria-label="Search variants">
    </div>
    <div class="table-wrap scroll"><table id="varTable"></table></div>
  </div>

  <div class="card" id="dataTables">
    <h2>Data tables</h2>
    <p class="sub">The CSV files saved with this page, shown here so you don't have to download them. They cover all samples.</p>
    <div class="filters" style="margin:0 0 10px">
      <div class="tabs" role="tablist" id="tblTabs"></div>
      <input type="search" id="tblSearch" placeholder="Search this table…" aria-label="Search data table">
      <span style="flex:1"></span>
      <button type="button" id="tblDownload">Download this CSV</button>
    </div>
    <p class="sub" id="tblInfo"></p>
    <div class="table-wrap scroll"><table id="dataTable"></table></div>
  </div>

  <footer id="foot"></footer>
</div>
<div id="tip" role="tooltip"></div>

<script type="application/json" id="data">__DATA__</script>
<script>
"use strict";
const DATA = JSON.parse(document.getElementById("data").textContent);
const T = DATA.thresholds;

// ---- colours (read from the CSS tokens, so light/dark both work) ----
function theme() {
  const cs = getComputedStyle(document.documentElement);
  const v = n => (cs.getPropertyValue(n) || "").trim();
  const dark = window.matchMedia && matchMedia("(prefers-color-scheme: dark)").matches;
  return {
    accent: v("--accent") || "#2a78d6", wash: v("--accent-wash") || "rgba(42,120,214,.1)",
    grid: v("--grid") || "#e8e7e3", text: v("--text") || "#0b0b0b", text2: v("--text-2") || "#52514e",
    muted: v("--muted") || "#7c7b75", surface: v("--surface") || "#fcfcfb",
    dropped: dark ? "#5d5c57" : "#bdbcb5", good: "#0ca30c", warning: "#fab219",
  };
}
let C = theme();

// ---- small helpers ----
const $ = id => document.getElementById(id);
const fmt = n => n == null ? "–" : Number(n).toLocaleString("en-US");
const fmt1 = n => n == null ? "–" : Number(n).toLocaleString("en-US", {maximumFractionDigits: 1});
const pct = n => n == null ? "–" : Number(n).toFixed(1) + "%";
const esc = s => String(s ?? "").replace(/[&<>"]/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;"}[c]));
const SVGNS = "http://www.w3.org/2000/svg";
function el(tag, attrs, parent) {
  const e = document.createElementNS(SVGNS, tag);
  for (const k in attrs || {}) e.setAttribute(k, attrs[k]);
  if (parent) parent.appendChild(e);
  return e;
}
function text(parent, x, y, str, opts) {
  const t = el("text", Object.assign({x, y, fill: C.text2, "font-size": 11}, opts || {}), parent);
  t.textContent = str;
  return t;
}
// Axis with round tick steps (1, 2 or 5 × 10^n) so labels read 0 / 10 / 20, not 0 / 13 / 25
function axis(dataMax, target) {
  const raw = Math.max(dataMax, 1e-9) / target;
  const p = Math.pow(10, Math.floor(Math.log10(raw))), f = raw / p;
  const step = (f <= 1 ? 1 : f <= 2 ? 2 : f <= 5 ? 5 : 10) * p;
  const max = Math.ceil(dataMax / step) * step || step;
  const t = []; for (let v = 0; v <= max + step / 2; v += step) t.push(+v.toFixed(10));
  return {max, ticks: t};
}
const compact = v => v >= 1e6 ? (v / 1e6).toFixed(v % 1e6 ? 1 : 0) + "M" : v >= 1e4 ? (v / 1e3).toFixed(0) + "K" : fmt(Math.round(v));
function width(id) { return Math.max(280, $(id).clientWidth || 900); }

// ---- tooltip ----
const tip = $("tip");
function showTip(html, evt) {
  tip.innerHTML = html; tip.style.display = "block";
  const pad = 14, r = tip.getBoundingClientRect();
  let x = evt.clientX + pad, y = evt.clientY + pad;
  if (x + r.width > innerWidth) x = evt.clientX - r.width - pad;
  if (y + r.height > innerHeight) y = evt.clientY - r.height - pad;
  tip.style.left = x + "px"; tip.style.top = y + "px";
}
function hideTip() { tip.style.display = "none"; }

// ---- state ----
let current = 0;
const S = () => DATA.samples[current];

// ===================== header, filters, samples table =====================
function renderMeta() {
  const r = DATA.reference;
  const len = r.contigs.reduce((a, c) => a + c.length, 0);
  $("meta").textContent = `${DATA.samples.length} sample${DATA.samples.length === 1 ? "" : "s"} · reference ${r.accession}` +
    (r.description ? ` (${r.description.replace(r.accession, "").trim()})` : "") +
    ` · ${fmt(len)} bp · generated ${DATA.generated}`;
  $("foot").innerHTML = `Run folder: <code>${esc(DATA.run)}</code><br>` +
    `CSV files for Power BI / Excel are saved next to this page in <code>07_dashboard/</code>: ` +
    `<code>sample_summary.csv</code>, <code>all_samples_variants.csv</code>, <code>all_samples_calls.csv</code>, ` +
    `<code>r_statistics.csv</code> and <code>variant_windows.csv</code>.`;
  const sel = $("sampleSel");
  sel.innerHTML = DATA.samples.map((s, i) => `<option value="${i}">${esc(s.name)}</option>`).join("");
  sel.onchange = () => select(+sel.value);
}

function renderSamplesTable() {
  const cols = [
    ["Sample", s => esc(s.name), ""],
    ["Reads after QC", s => fmt(s.qc?.reads_pass), "num"],
    ["Mapped", s => pct(s.align?.mapped_pct), "num"],
    ["Mean depth", s => fmt1(s.cov?.meandepth) + "×", "num"],
    ["Genome covered", s => pct(s.cov?.covered_pct), "num"],
    ["Raw calls", s => fmt(s.calls.length), "num"],
    ["Trusted variants", s => fmt(s.variants.length), "num"],
    ["Status", s => statusBadge(checks(s)), ""],
  ];
  $("samplesTable").innerHTML =
    "<thead><tr>" + cols.map(c => `<th class="${c[2]}">${c[0]}</th>`).join("") + "</tr></thead><tbody>" +
    DATA.samples.map((s, i) => `<tr class="clickable${i === current ? " selected" : ""}" data-i="${i}">` +
      cols.map(c => `<td class="${c[2]}">${c[1](s)}</td>`).join("") + "</tr>").join("") + "</tbody>";
  $("samplesTable").querySelectorAll("tr.clickable").forEach(tr => tr.onclick = () => select(+tr.dataset.i));
}

// ===================== quality checks =====================
// WHY these four: each one, if it fails, makes the variant list less reliable.
function checks(s) {
  const out = [];
  const qc = s.qc, al = s.align, cv = s.cov;
  if (qc) {
    const kept = 100 * qc.reads_pass / qc.reads_raw;
    out.push({ok: kept >= 80, name: "Reads passing QC", value: pct(kept),
      why: kept >= 80 ? "Most reads were good quality." : "Many reads were removed for low quality – check the fastp report."});
  }
  if (al) out.push({ok: al.mapped_pct >= 50, name: "Reads mapped to the reference", value: pct(al.mapped_pct),
      why: al.mapped_pct >= 50 ? "Most reads come from this genome." :
        "Most reads don't match the reference – often host or other organisms in the sample. Expected for clinical RNA-seq; worrying for a pure culture."});
  if (cv) {
    out.push({ok: cv.meandepth >= T.min_dp, name: "Mean read depth", value: fmt1(cv.meandepth) + "×",
      why: cv.meandepth >= T.min_dp ? `Above the ${T.min_dp}× needed to keep a variant.` :
        `Below ${T.min_dp}×, so many positions can't pass the depth filter.`});
    out.push({ok: cv.covered_pct >= 90, name: "Genome covered by at least one good-quality read", value: pct(cv.covered_pct),
      why: cv.covered_pct >= 90 ? "Almost the whole genome was seen." : "Parts of the genome have no reads – variants there can't be detected."});
  }
  // WHY: in a haploid genome real variants are in (nearly) all reads. Several
  //      calls at 20–80% point to a mixed sample, contamination or mapping errors.
  if (s.calls.some(c => c.af != null)) {
    const mixed = s.calls.filter(c => c.af != null && c.af >= 0.2 && c.af < T.min_af).length;
    out.push({ok: mixed <= 5, name: `Mixed calls (AF 20–${Math.round(T.min_af * 100)}%)`, value: fmt(mixed),
      why: mixed <= 5 ? "Variants are in (nearly) all reads, as expected for a single haploid genome." :
        "Many positions have a mix of reads – possible mixed infection, contamination, or reads from a related genome. These were flagged LowAF."});
  }
  return out;
}
function statusBadge(list) {
  const bad = list.filter(c => !c.ok).length;
  return bad ? `<span class="badge"><span class="dot" style="background:${C.warning}"></span>${bad} warning${bad > 1 ? "s" : ""}</span>`
             : `<span class="badge"><span class="dot" style="background:${C.good}"></span>All checks passed</span>`;
}
function renderChecks() {
  const list = checks(S());
  const ann = S().annotated
    ? `<li><span class="badge"><span class="dot" style="background:${C.good}"></span>✓ Done</span><div><b>Annotation</b><div class="why">Genes and effects were added with SnpEff.</div></div></li>`
    : `<li><span class="badge"><span class="dot" style="background:${C.muted}"></span>– Skipped</span><div><b>Annotation</b><div class="why">No SnpEff database is set up for ${esc(DATA.reference.accession)} (SNPEFF_DB), so gene and effect columns are empty.</div></div></li>`;
  $("checks").innerHTML = list.map(c =>
    `<li><span class="badge"><span class="dot" style="background:${c.ok ? C.good : C.warning}"></span>${c.ok ? "✓ OK" : "! Check"}</span>` +
    `<div><b>${c.name}: ${c.value}</b><div class="why">${c.why}</div></div></li>`).join("") + ann;
}

// ===================== tiles =====================
function renderTiles() {
  const s = S(), snps = s.variants.filter(v => v.REF.length === 1 && v.ALT.length === 1).length;
  const tiles = [
    ["Reads after QC", fmt(s.qc?.reads_pass), s.qc ? `of ${fmt(s.qc.reads_raw)} sequenced` : ""],
    ["Mapped to reference", pct(s.align?.mapped_pct), s.align ? `${fmt(s.align.mapped)} reads` : ""],
    ["Mean read depth", fmt1(s.cov?.meandepth) + "×", `good-quality reads · min. ${T.min_dp}× to keep a variant`],
    ["Genome covered", pct(s.cov?.covered_pct), "by at least one good-quality read"],
    ["Trusted variants", fmt(s.variants.length), `of ${fmt(s.calls.length)} calls · ${fmt(snps)} SNPs · ${fmt(s.variants.length - snps)} other`],
  ];
  $("tiles").innerHTML = tiles.map(t =>
    `<div class="tile"><div class="label">${t[0]}</div><div class="value">${t[1]}</div><div class="note">${t[2]}</div></div>`).join("");
}

// ===================== funnel (horizontal bars) =====================
function renderFunnel() {
  const s = S(), calls = s.calls;
  // Cumulative: each bar applies one more filter (AF is skipped for old runs without it)
  const okQ = c => c.qual >= T.min_qual, okD = c => okQ(c) && c.dp >= T.min_dp;
  const steps = [
    ["Raw calls (bcftools)", calls.length],
    [`+ QUAL ≥ ${T.min_qual}`, calls.filter(okQ).length],
    [`+ depth ≥ ${T.min_dp} good-quality reads`, calls.filter(okD).length],
  ];
  if (calls.some(c => c.af != null))
    steps.push([`+ AF ≥ ${Math.round(T.min_af * 100)}% of reads (= trusted)`, calls.filter(c => c.kept).length]);
  $("funnelSub").textContent = `Each filter removes calls that are more likely to be errors. Removed calls are flagged, not deleted (see "All calls" table).`;
  const W = width("funnel"), bar = 22, gap = 22, labelH = 16, H = steps.length * (bar + gap + labelH);
  const svg = el("svg", {width: W, height: H, role: "img", "aria-label": "Filter funnel"});
  const max = Math.max(1, steps[0][1]), right = 56, plotW = W - right;
  steps.forEach(([label, n], i) => {
    const y = i * (bar + gap + labelH);
    text(svg, 0, y + 12, label, {fill: C.text2, "font-size": 12});
    const w = Math.max(4, plotW * n / max), by = y + labelH;
    el("rect", {x: 0, y: by, width: plotW, height: bar, rx: 4, fill: C.grid}, svg);
    // bar with a rounded data end, square at the baseline
    el("path", {d: `M0,${by} H${w - 4} a4,4 0 0 1 4,4 V${by + bar - 4} a4,4 0 0 1 -4,4 H0 Z`, fill: C.accent}, svg);
    text(svg, Math.min(w + 8, plotW + 6), by + bar / 2 + 4, fmt(n), {fill: C.text, "font-size": 12, "font-weight": 600});
    const hit = el("rect", {x: 0, y: by - 4, width: W, height: bar + 8, fill: "transparent"}, svg);
    hit.onmousemove = e => showTip(`<b>${esc(label)}</b><br>${fmt(n)} calls (${(100 * n / max).toFixed(0)}% of raw)`, e);
    hit.onmouseleave = hideTip;
  });
  $("funnel").replaceChildren(svg);
}

// ===================== coverage along genome + variant ticks =====================
function renderCoverage() {
  const s = S(), prof = s.profile, box = $("coverage");
  if (!prof) { box.innerHTML = '<p class="empty">No alignment found.</p>'; return; }
  const glen = DATA.reference.contigs.reduce((a, c) => a + c.length, 0);
  const st = s.settings || {};
  $("covSub").textContent = (st.min_mapq == null ? "Average reads " :
    `Average good-quality reads (mapping quality ≥ ${st.min_mapq}, base quality ≥ ${st.min_baseq}) `) +
    `covering each part of the genome (${fmt(prof.bin_size)} bp per point). ` +
    `Blue ticks below are the trusted variants. Below the ${T.min_dp}× line, variants can't pass the depth filter.`;
  const W = width("coverage"), H = 260, m = {l: 44, r: 12, t: 10, b: 58};
  const pw = W - m.l - m.r, ph = H - m.t - m.b;
  const ya = axis(Math.max(T.min_dp * 1.2, ...prof.depth), 4), ymax = ya.max;
  const x = g => m.l + pw * g / glen, y = d => m.t + ph - ph * d / ymax;
  const svg = el("svg", {width: W, height: H, role: "img", "aria-label": "Read depth along the genome"});
  ya.ticks.forEach(t => {
    el("line", {x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), stroke: C.grid, "stroke-width": 1}, svg);
    text(svg, m.l - 6, y(t) + 4, compact(t) + "×", {"text-anchor": "end"});
  });
  axis(glen, 5).ticks.filter(t => t <= glen).forEach(t =>
    text(svg, x(t), m.t + ph + 16, compact(t), {"text-anchor": t === 0 ? "start" : "middle"}));
  text(svg, m.l + pw / 2, H - 4, "Position in genome (bp)", {"text-anchor": "middle"});
  const pts = prof.depth.map((d, i) => [x((i + 0.5) * prof.bin_size), y(d)]);
  const line = pts.map((p, i) => (i ? "L" : "M") + p[0].toFixed(1) + "," + p[1].toFixed(1)).join("");
  el("path", {d: line + `L${pts[pts.length - 1][0]},${y(0)}L${pts[0][0]},${y(0)}Z`, fill: C.wash}, svg);
  el("path", {d: line, fill: "none", stroke: C.accent, "stroke-width": 2, "stroke-linejoin": "round", "stroke-linecap": "round"}, svg);
  // minimum-depth reference line
  el("line", {x1: m.l, x2: W - m.r, y1: y(T.min_dp), y2: y(T.min_dp), stroke: C.muted, "stroke-width": 1}, svg);
  // label on a small surface-coloured box so the data line doesn't run through the text
  const lbl = `min. depth ${T.min_dp}×`, lw = lbl.length * 6.4 + 8;
  el("rect", {x: W - m.r - lw, y: y(T.min_dp) - 18, width: lw, height: 15, rx: 3, fill: C.surface, opacity: 0.92}, svg);
  text(svg, W - m.r - 4, y(T.min_dp) - 7, lbl, {"text-anchor": "end", fill: C.text2});
  // variant tick strip
  const ty = m.t + ph + 24;
  el("line", {x1: m.l, x2: W - m.r, y1: ty + 5, y2: ty + 5, stroke: C.grid}, svg);
  s.calls.filter(c => c.kept).forEach(c => {
    el("line", {x1: x(c.gpos), x2: x(c.gpos), y1: ty, y2: ty + 10, stroke: C.accent, "stroke-width": 1.5}, svg);
  });
  // crosshair + tooltip
  const cross = el("line", {y1: m.t, y2: m.t + ph, stroke: C.muted, "stroke-width": 1, visibility: "hidden"}, svg);
  const dot = el("circle", {r: 4, fill: C.accent, stroke: C.surface, "stroke-width": 2, visibility: "hidden"}, svg);
  const hit = el("rect", {x: m.l, y: m.t, width: pw, height: ph + 40, fill: "transparent"}, svg);
  hit.onmousemove = e => {
    const r = svg.getBoundingClientRect(), gx = (e.clientX - r.left - m.l) / pw * glen;
    const i = Math.max(0, Math.min(prof.depth.length - 1, Math.floor(gx / prof.bin_size)));
    const px = x((i + 0.5) * prof.bin_size);
    cross.setAttribute("x1", px); cross.setAttribute("x2", px); cross.setAttribute("visibility", "visible");
    dot.setAttribute("cx", px); dot.setAttribute("cy", y(prof.depth[i])); dot.setAttribute("visibility", "visible");
    const lo = i * prof.bin_size + 1, hi = Math.min(glen, (i + 1) * prof.bin_size);
    const nv = s.calls.filter(c => c.kept && c.gpos >= lo && c.gpos <= hi).length;
    showTip(`<b>${fmt(lo)}–${fmt(hi)} bp</b><br>Mean depth: ${fmt1(prof.depth[i])}×<br>Trusted variants here: ${nv}`, e);
  };
  hit.onmouseleave = () => { cross.setAttribute("visibility", "hidden"); dot.setAttribute("visibility", "hidden"); hideTip(); };
  box.replaceChildren(svg);
}

// ===================== QUAL vs depth scatter =====================
function renderScatter() {
  const s = S(), calls = s.calls, box = $("scatter");
  $("lgKept").style.background = C.accent; $("lgDrop").style.background = C.dropped;
  $("scatterSub").textContent = `Each dot is one raw call. Kept = QUAL ≥ ${T.min_qual}, depth ≥ ${T.min_dp} (top-right of the lines) ` +
    `and AF ≥ ${Math.round(T.min_af * 100)}%. Hover a dot for its filter reason.`;
  if (!calls.length) { box.innerHTML = '<p class="empty">No variant calls.</p>'; return; }
  const W = width("scatter"), H = 300, m = {l: 44, r: 12, t: 10, b: 40};
  const pw = W - m.l - m.r, ph = H - m.t - m.b;
  const xa = axis(Math.max(T.min_dp * 1.5, ...calls.map(c => c.dp)), 5), xmax = xa.max;
  const ya = axis(Math.max(T.min_qual * 1.5, ...calls.map(c => c.qual)), 4), ymax = ya.max;
  const x = v => m.l + pw * v / xmax, y = v => m.t + ph - ph * v / ymax;
  const svg = el("svg", {width: W, height: H, role: "img", "aria-label": "QUAL against read depth"});
  ya.ticks.forEach(t => {
    el("line", {x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), stroke: C.grid}, svg);
    text(svg, m.l - 6, y(t) + 4, fmt(t), {"text-anchor": "end"});
  });
  xa.ticks.forEach(t => text(svg, x(t), m.t + ph + 16, fmt(t), {"text-anchor": "middle"}));
  text(svg, m.l + pw / 2, H - 4, "Read depth (DP)", {"text-anchor": "middle"});
  text(svg, 12, m.t + ph / 2, "QUAL", {"text-anchor": "middle", transform: `rotate(-90 12 ${m.t + ph / 2})`});
  el("line", {x1: x(T.min_dp), x2: x(T.min_dp), y1: m.t, y2: m.t + ph, stroke: C.muted}, svg);
  el("line", {x1: m.l, x2: W - m.r, y1: y(T.min_qual), y2: y(T.min_qual), stroke: C.muted}, svg);
  // filtered first, kept on top
  [...calls].sort((a, b) => a.kept - b.kept).forEach(c => {
    const d = el("circle", {cx: x(c.dp), cy: y(c.qual), r: 4, fill: c.kept ? C.accent : C.dropped,
      stroke: C.surface, "stroke-width": 1.5}, svg);
    const hit = el("circle", {cx: x(c.dp), cy: y(c.qual), r: 8, fill: "transparent"}, svg);
    hit.onmousemove = e => showTip(`<b>${esc(c.chrom)}:${fmt(c.pos)}</b> ${esc(c.ref)}→${esc(c.alt)}<br>` +
      `QUAL ${fmt1(c.qual)} · depth ${c.dp}` + (c.af != null ? ` · AF ${Math.round(c.af * 100)}%` : "") +
      `<br>${c.kept ? "Kept (PASS)" : "Filtered out: " + esc(c.filter.replace(/;/g, ", "))}`, e);
    hit.onmouseleave = hideTip;
  });
  box.replaceChildren(svg);
}

// ===================== change types / effects (columns) =====================
function renderTypes() {
  const s = S(), box = $("types");
  let cats;
  if (s.annotated) {
    $("typesTitle").textContent = "Effect on genes";
    $("typesSub").textContent = "What the trusted variants do, according to SnpEff (most severe effect per variant).";
    const counts = {};
    s.variants.forEach(v => { const k = (v.Effect || "unknown").split("&")[0].replace(/_/g, " "); counts[k] = (counts[k] || 0) + 1; });
    cats = Object.entries(counts).sort((a, b) => b[1] - a[1]).slice(0, 10);
  } else {
    $("typesTitle").textContent = "Types of change";
    const order = ["A>G", "G>A", "C>T", "T>C", "A>C", "C>A", "A>T", "T>A", "G>C", "C>G", "G>T", "T>G"];
    const counts = Object.fromEntries(order.map(k => [k, 0])); let indel = 0, ts = 0, tv = 0;
    s.variants.forEach(v => {
      if (v.REF.length === 1 && v.ALT.length === 1) {
        const k = v.REF + ">" + v.ALT; if (k in counts) counts[k]++;
        (order.indexOf(k) < 4 ? ts++ : tv++);
      } else indel++;
    });
    cats = order.map(k => [k, counts[k]]).concat([["Indel/other", indel]]);
    $("typesSub").textContent = `Base changes among the trusted variants. The first four are transitions (Ts = ${ts}), the rest transversions (Tv = ${tv})` +
      (tv ? `; Ts/Tv = ${(ts / tv).toFixed(2)}.` : ".");
  }
  if (!cats.length || !s.variants.length) { box.innerHTML = '<p class="empty">No trusted variants.</p>'; return; }
  const W = width("types"), H = 300, m = {l: 40, r: 8, t: 18, b: s.annotated ? 90 : 40};
  const pw = W - m.l - m.r, ph = H - m.t - m.b, band = pw / cats.length, bw = Math.min(24, band * 0.6);
  const ya = axis(Math.max(...cats.map(c => c[1])), 4), ymax = ya.max;
  const y = v => m.t + ph - ph * v / ymax;
  const svg = el("svg", {width: W, height: H, role: "img", "aria-label": s.annotated ? "Effects" : "Base change types"});
  ya.ticks.forEach(t => {
    el("line", {x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), stroke: C.grid}, svg);
    text(svg, m.l - 6, y(t) + 4, fmt(t), {"text-anchor": "end"});
  });
  cats.forEach(([k, n], i) => {
    const cx = m.l + band * (i + 0.5), x0 = cx - bw / 2, top = y(n), base = y(0);
    if (n > 0) el("path", {d: `M${x0},${base} V${top + 4} a4,4 0 0 1 4,-4 H${x0 + bw - 4} a4,4 0 0 1 4,4 V${base} Z`, fill: C.accent}, svg);
    if (n > 0) text(svg, cx, top - 4, fmt(n), {"text-anchor": "middle", fill: C.text, "font-size": 11});
    const lbl = text(svg, cx, base + 14, k.length > 18 ? k.slice(0, 17) + "…" : k,
      s.annotated ? {"text-anchor": "end", transform: `rotate(-35 ${cx} ${base + 14})`} : {"text-anchor": "middle", "font-size": 10});
    const hit = el("rect", {x: cx - band / 2, y: m.t, width: band, height: ph, fill: "transparent"}, svg);
    hit.onmousemove = e => showTip(`<b>${esc(k)}</b><br>${fmt(n)} variant${n === 1 ? "" : "s"}`, e);
    hit.onmouseleave = hideTip;
  });
  box.replaceChildren(svg);
}

// ===================== variant table (search + sort) =====================
let sortKey = "POS", sortDir = 1;
function renderVariants() {
  const s = S(), q = $("varSearch").value.trim().toLowerCase();
  const cols = [["POS", "Position", "num"], ["REF", "Ref", ""], ["ALT", "Alt", ""], ["QUAL", "QUAL", "num"],
                ["DP", "Depth", "num"], ["AF", "AF", "num"], ["Effect", "Effect", ""], ["Impact", "Impact", ""], ["Gene", "Gene", ""], ["HGVS_p", "Protein change", ""]];
  let rows = s.variants.filter(v => !q || Object.values(v).join(" ").toLowerCase().includes(q));
  const numeric = ["POS", "QUAL", "DP", "AF"].includes(sortKey);
  rows.sort((a, b) => (numeric ? (+a[sortKey] - +b[sortKey]) : String(a[sortKey] || "").localeCompare(String(b[sortKey] || ""))) * sortDir);
  $("varSub").textContent = `${fmt(rows.length)} of ${fmt(s.variants.length)} variants shown. Click a column to sort.` +
    (s.annotated ? "" : " Effect/Impact/Gene are empty because annotation was skipped.");
  $("varTable").innerHTML = "<thead><tr>" + cols.map(c =>
      `<th class="${c[2]}" data-k="${c[0]}">${c[1]}${sortKey === c[0] ? (sortDir > 0 ? " ▲" : " ▼") : ""}</th>`).join("") +
    "</tr></thead><tbody>" + rows.map(v => "<tr>" + cols.map(c => {
      let val = v[c[0]] ?? "";
      if (c[0] === "POS") val = fmt(val);
      if (c[0] === "QUAL") val = fmt1(val);
      if (c[0] === "AF" && val !== "") val = Math.round(+val * 100) + "%";
      return `<td class="${c[2]}">${esc(val)}</td>`;
    }).join("") + "</tr>").join("") + "</tbody>";
  $("varTable").querySelectorAll("th").forEach(th => th.onclick = () => {
    const k = th.dataset.k; sortDir = sortKey === k ? -sortDir : 1; sortKey = k; renderVariants();
  });
}
$("varSearch").oninput = renderVariants;

// ===================== Data tables: the CSV files, viewable in the page =====================
// WHY: opening a CSV from a browser usually just downloads it. These are the
//      exact same rows (built once in Python), shown as searchable tables.
const TABLES = DATA.tables || {};
const ROW_LIMIT = 1000;   // WHY: keeps the page fast; search narrows large tables
let tbl = "summary", tblSort = null, tblDir = 1;
function isNumCol(rows, i) {
  return rows.length > 0 && rows.every(r => r[i] === null || r[i] === "" || (typeof r[i] !== "boolean" && !isNaN(+r[i])));
}
const cellText = v => v === true ? "True" : v === false ? "False" : v;   // same as in the CSV
const cellNum = v => Number(v).toLocaleString("en-US", {maximumFractionDigits: 3});   // keep the CSV's decimals
function renderDataTable() {
  const t = TABLES[tbl]; if (!t) return;
  document.querySelectorAll(".tabs button").forEach(b => b.setAttribute("aria-selected", b.dataset.t === tbl));
  const q = $("tblSearch").value.trim().toLowerCase();
  const num = t.columns.map((_, i) => isNumCol(t.rows, i));
  let rows = t.rows.filter(r => !q || r.join(" ").toLowerCase().includes(q));
  if (tblSort !== null) {
    const i = tblSort;
    rows = [...rows].sort((a, b) => (num[i] ? (+a[i] || 0) - (+b[i] || 0) : String(a[i] ?? "").localeCompare(String(b[i] ?? ""))) * tblDir);
  }
  const shown = rows.slice(0, ROW_LIMIT);
  $("tblInfo").textContent = `${t.file} · ${fmt(rows.length)} of ${fmt(t.rows.length)} rows` +
    (rows.length > ROW_LIMIT ? ` (first ${fmt(ROW_LIMIT)} shown – search to narrow down)` : "") + " · click a column to sort";
  $("dataTable").innerHTML = "<thead><tr>" + t.columns.map((c, i) =>
      `<th class="${num[i] ? "num" : ""}" data-i="${i}">${esc(c.replace(/_/g, " "))}${tblSort === i ? (tblDir > 0 ? " ▲" : " ▼") : ""}</th>`).join("") +
    "</tr></thead><tbody>" + shown.map(r => "<tr>" + r.map((v, i) =>
      `<td class="${num[i] ? "num" : ""}">${esc(v === null ? "" : (num[i] && v !== "" ? cellNum(v) : cellText(v)))}</td>`).join("") + "</tr>").join("") + "</tbody>";
  $("dataTable").querySelectorAll("th").forEach(th => th.onclick = () => {
    const i = +th.dataset.i; tblDir = tblSort === i ? -tblDir : 1; tblSort = i; renderDataTable();
  });
}
$("tblTabs").innerHTML = Object.entries(TABLES).map(([k, t]) =>
  `<button type="button" role="tab" data-t="${k}">${esc(t.label || t.file)}</button>`).join("");
document.querySelectorAll(".tabs button").forEach(b => b.onclick = () => { tbl = b.dataset.t; tblSort = null; renderDataTable(); });
$("tblSearch").oninput = renderDataTable;
// Download is still there for when you do want the file (whole table, not just the search result)
$("tblDownload").onclick = () => {
  const t = TABLES[tbl]; if (!t) return;
  const cell = v => { const s = String(v ?? ""); return /[",\n]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s; };
  const csv = [t.columns, ...t.rows].map(r => r.map(v => cell(cellText(v))).join(",")).join("\n");
  const a = document.createElement("a");
  a.href = URL.createObjectURL(new Blob([csv], {type: "text/csv"})); a.download = t.file; a.click();
  setTimeout(() => URL.revokeObjectURL(a.href), 1000);
};

// ===================== wiring =====================
// ===================== AF histogram =====================
// WHY: a clean haploid sample has almost all calls near 100%. A second bump
//      in the middle means mixed reads (mixed infection / contamination).
function renderAF() {
  const s = S(), box = $("afhist"), calls = s.calls.filter(c => c.af != null);
  $("afSub").textContent = `Share of good-quality reads showing the variant, for every call. ` +
    `Calls left of the ${Math.round(T.min_af * 100)}% line are flagged LowAF.`;
  if (!calls.length) { box.innerHTML = '<p class="empty">No allele frequencies in this run (made before AF was recorded – re-run to see this).</p>'; return; }
  const bins = Array.from({length: 10}, (_, i) => ({lo: i / 10, n: 0}));
  calls.forEach(c => bins[Math.min(9, Math.floor(c.af * 10))].n++);
  const W = width("afhist"), H = 260, m = {l: 40, r: 8, t: 18, b: 40};
  const pw = W - m.l - m.r, ph = H - m.t - m.b, band = pw / 10, bw = Math.min(24, band * 0.7);
  const ya = axis(Math.max(...bins.map(b => b.n)), 4), y = v => m.t + ph - ph * v / ya.max, x = f => m.l + pw * f;
  const svg = el("svg", {width: W, height: H, role: "img", "aria-label": "Allele frequency histogram"});
  ya.ticks.forEach(t => { el("line", {x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), stroke: C.grid}, svg);
    text(svg, m.l - 6, y(t) + 4, fmt(t), {"text-anchor": "end"}); });
  [0, 0.2, 0.4, 0.6, 0.8, 1].forEach(f => text(svg, x(f), m.t + ph + 16, Math.round(f * 100) + "%",
    {"text-anchor": f === 0 ? "start" : f === 1 ? "end" : "middle"}));   // end labels stay inside the chart
  text(svg, m.l + pw / 2, H - 4, "Allele frequency (AF)", {"text-anchor": "middle"});
  bins.forEach((b, i) => {
    const cx = x(b.lo + 0.05), x0 = cx - bw / 2, top = y(b.n), base = y(0);
    const fill = b.lo + 0.1 <= T.min_af ? C.dropped : C.accent;
    if (b.n > 0) el("path", {d: `M${x0},${base} V${top + 4} a4,4 0 0 1 4,-4 H${x0 + bw - 4} a4,4 0 0 1 4,4 V${base} Z`, fill}, svg);
    if (b.n > 0) text(svg, cx, top - 4, fmt(b.n), {"text-anchor": "middle", fill: C.text});
    const hit = el("rect", {x: cx - band / 2, y: m.t, width: band, height: ph, fill: "transparent"}, svg);
    hit.onmousemove = e => showTip(`<b>AF ${Math.round(b.lo * 100)}–${Math.round(b.lo * 100) + 10}%</b><br>${fmt(b.n)} call${b.n === 1 ? "" : "s"}`, e);
    hit.onmouseleave = hideTip;
  });
  el("line", {x1: x(T.min_af), x2: x(T.min_af), y1: m.t, y2: m.t + ph, stroke: C.muted}, svg);
  const lbl = `min. AF ${Math.round(T.min_af * 100)}%`, lw = lbl.length * 6.4 + 8;
  el("rect", {x: x(T.min_af) - lw - 4, y: m.t, width: lw, height: 15, rx: 3, fill: C.surface, opacity: 0.92}, svg);
  text(svg, x(T.min_af) - 8, m.t + 11, lbl, {"text-anchor": "end"});
  box.replaceChildren(svg);
}

// ===================== filter reasons =====================
// WHY: shows which rule removed most calls – e.g. mostly LowDepth means
//      "sequence deeper", mostly LowAF means "check for a mixed sample".
function renderReasons() {
  const s = S(), box = $("reasons"), calls = s.calls;
  const rows = [
    ["LowQual", `QUAL below ${T.min_qual}`],
    ["LowDepth", `fewer than ${T.min_dp} good-quality reads`],
    ["LowAF", `fewer than ${Math.round(T.min_af * 100)}% of reads show it`],
  ].map(([k, d]) => [k, d, calls.filter(c => c.filter.split(";").includes(k)).length]);
  const flagged = calls.filter(c => !c.kept).length;
  $("reasonSub").textContent = `${fmt(flagged)} of ${fmt(calls.length)} calls were flagged. A call can fail more than one rule.`;
  const W = width("reasons"), bar = 22, gap = 22, labelH = 16, H = rows.length * (bar + gap + labelH);
  const svg = el("svg", {width: W, height: H, role: "img", "aria-label": "Filter reasons"});
  const max = Math.max(1, calls.length), plotW = W - 56;
  rows.forEach(([k, d, n], i) => {
    const y0 = i * (bar + gap + labelH), by = y0 + labelH, w = n ? Math.max(4, plotW * n / max) : 0;
    text(svg, 0, y0 + 12, `${k} – ${d}`, {"font-size": 12});
    el("rect", {x: 0, y: by, width: plotW, height: bar, rx: 4, fill: C.grid}, svg);
    if (n) el("path", {d: `M0,${by} H${w - 4} a4,4 0 0 1 4,4 V${by + bar - 4} a4,4 0 0 1 -4,4 H0 Z`, fill: C.dropped}, svg);
    text(svg, Math.min(w + 8, plotW + 6), by + bar / 2 + 4, fmt(n), {fill: C.text, "font-size": 12, "font-weight": 600});
    const hit = el("rect", {x: 0, y: by - 4, width: W, height: bar + 8, fill: "transparent"}, svg);
    hit.onmousemove = e => showTip(`<b>${k}</b><br>${fmt(n)} of ${fmt(calls.length)} calls (${(100 * n / max).toFixed(0)}%)`, e);
    hit.onmouseleave = hideTip;
  });
  box.replaceChildren(svg);
}

// ===================== calling settings =====================
function renderSettings() {
  const st = S().settings || {};
  // ok: true = ✓ On, false = ! Off, null = unknown (run made before settings were saved)
  const item = (ok, name, why) => {
    const [col, lab] = ok == null ? [C.muted, "? Unknown"] : ok ? [C.good, "✓ On"] : [C.warning, "! Off"];
    return `<li><span class="badge"><span class="dot" style="background:${col}"></span>${lab}</span>` +
      `<div><b>${name}</b><div class="why">${ok == null ? "Not recorded for this run – re-run to record it. " : ""}${why}</div></div></li>`;
  };
  const v = (x, test) => x == null ? null : test(x);
  $("settings").innerHTML = [
    item(v(st.min_mapq, x => x >= 20), `Mapping quality ≥ ${st.min_mapq ?? "?"}`, "Ignores reads that could belong in several places."),
    item(v(st.min_baseq, x => x >= 20), `Base quality ≥ ${st.min_baseq ?? "?"}`, "Ignores blurry bases (Q20 = 99% accurate)."),
    item(v(st.max_depth, x => x >= 1000), `Up to ${st.max_depth == null ? "?" : fmt(st.max_depth)} reads per position`, "Deep data isn't silently cut (the old default was 250)."),
    item(st.haploid, "Haploid calling (--ploidy 1)", "One genome copy, as in bacteria and viruses."),
    item(st.normalized, "Indels normalized (bcftools norm)", "Each variant is written one standard way, so samples compare correctly."),
    item(S().calls.some(c => c.af != null), `Filters: QUAL ≥ ${T.min_qual} · depth ≥ ${T.min_dp} · AF ≥ ${Math.round(T.min_af * 100)}%`,
      "Failing calls are flagged with the reason (FILTER column), not deleted."),
  ].join("");
}

// ===================== statistics from R (step 9b) =====================
// WHY: counts say how many variants there are; the tests say whether the
//      filter really separates good from bad calls, whether the base changes
//      look like real mutations (Ts/Tv), and where variants cluster.
function renderR() {
  const r = S().r, card = $("rCard");
  card.hidden = !r;            // runs made before the R step: hide the section
  if (!r) return;
  $("rTable").innerHTML = "<thead><tr><th>Test</th><th>Result</th><th>What it means</th></tr></thead><tbody>" +
    r.stats.map(x => `<tr><td>${esc(x.Metric)}</td><td class="num">${esc(x.Value ?? "")}</td><td class="why">${esc(x.Interpretation ?? "")}</td></tr>`).join("") +
    "</tbody>";
  renderWindows();
}
// Bar per genome window: trusted variants found vs the count expected if they
// were spread evenly (dashed line). Hotspots (Poisson test, BH-corrected) in accent colour.
function renderWindows() {
  const r = S().r, box = $("windows"), w = (r && r.windows) || [];
  $("lgHot").style.background = C.accent; $("lgWin").style.background = C.dropped; $("lgExp").style.background = C.muted;
  if (!w.length) { box.innerHTML = '<p class="empty">No window table.</p>'; return; }
  const hot = w.filter(x => x.Hotspot).length, size = w[0].Window_end - w[0].Window_start + 1;
  $("winSub").textContent = `Trusted variants per ${fmt(size)} bp window. ${hot} window${hot === 1 ? "" : "s"} ha${hot === 1 ? "s" : "ve"} more than chance would give.`;
  const W = width("windows"), H = 240, m = {l: 36, r: 8, t: 14, b: 36};
  const pw = W - m.l - m.r, ph = H - m.t - m.b, band = pw / w.length, bw = Math.max(2, band * 0.75);
  const ya = axis(Math.max(...w.map(x => Math.max(x.Variants, x.Expected))), 4), y = v => m.t + ph - ph * v / ya.max;
  const svg = el("svg", {width: W, height: H, role: "img", "aria-label": "Variants per genome window"});
  ya.ticks.forEach(t => { el("line", {x1: m.l, x2: W - m.r, y1: y(t), y2: y(t), stroke: C.grid}, svg);
    text(svg, m.l - 6, y(t) + 4, fmt(t), {"text-anchor": "end"}); });
  const genome = w[w.length - 1].Window_end;
  [0, 0.5, 1].forEach(f => text(svg, m.l + pw * f, m.t + ph + 16, fmt(Math.round(genome * f)),
    {"text-anchor": f === 0 ? "start" : f === 1 ? "end" : "middle"}));
  text(svg, m.l + pw / 2, H - 4, "Genome position (bp)", {"text-anchor": "middle"});
  w.forEach((x, i) => {
    const cx = m.l + band * (i + 0.5), top = y(x.Variants), base = y(0);
    if (x.Variants > 0) el("rect", {x: cx - bw / 2, y: top, width: bw, height: base - top, rx: Math.min(3, bw / 2),
      fill: x.Hotspot ? C.accent : C.dropped}, svg);
    const hit = el("rect", {x: cx - band / 2, y: m.t, width: band, height: ph, fill: "transparent"}, svg);
    hit.onmousemove = e => showTip(`<b>${fmt(x.Window_start)}–${fmt(x.Window_end)} bp</b><br>${x.Variants} variants (expected ${x.Expected})` +
      `<br>q = ${x.Q_value}${x.Hotspot ? " · <b>hotspot</b>" : ""}`, e);
    hit.onmouseleave = hideTip;
  });
  const ey = y(w[0].Expected);
  el("line", {x1: m.l, x2: W - m.r, y1: ey, y2: ey, stroke: C.muted, "stroke-dasharray": "4 3"}, svg);
  box.replaceChildren(svg);
}

function renderCharts() { renderFunnel(); renderCoverage(); renderScatter(); renderTypes(); renderAF(); renderReasons(); renderWindows(); }
function select(i) {
  current = i; $("sampleSel").value = i;
  renderSamplesTable(); renderTiles(); renderChecks(); renderCharts(); renderVariants(); renderSettings(); renderR();
}
if (!DATA.samples.length) {
  document.querySelector(".wrap").insertAdjacentHTML("beforeend", '<p class="empty">No finished samples found (no files in 06_reports/).</p>');
} else {
  renderMeta();
  if (DATA.samples.length < 2) $("samplesCard").hidden = true;
  select(0);
  renderDataTable();
  let t; addEventListener("resize", () => { clearTimeout(t); t = setTimeout(renderCharts, 150); });
  if (window.matchMedia) matchMedia("(prefers-color-scheme: dark)").addEventListener?.("change", () => { C = theme(); select(current); });
}
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Command line
# ---------------------------------------------------------------------------
if __name__ == "__main__":
    p = argparse.ArgumentParser(description="Build the results dashboard for a run folder.")
    p.add_argument("run_dir", help="output folder of a run (contains 06_reports/)")
    p.add_argument("--min-qual", type=float, default=30.0)
    p.add_argument("--min-dp", type=int, default=15)
    p.add_argument("--min-af", type=float, default=0.8)
    a = p.parse_args()
    out = build_dashboard(Path(a.run_dir), a.min_qual, a.min_dp, a.min_af)
    print(f"Dashboard written → {out}")
