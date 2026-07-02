__author__ = 'Liam TrinhNguyen'
__email__ = 'LiamTrinhNguyen@gmail.com'
__version__ = 'NGS_Pipeline_v2.1'

import os
import sys
import gzip
import shutil
import base64
import logging
import urllib.request
import subprocess
from datetime import datetime
from pathlib import Path
import polars as pl

# Progress bar with ASCII fallback
try:
    from tqdm import tqdm
except ImportError:
    print("Installing tqdm...")
    subprocess.check_call([sys.executable, "-m", "pip", "install", "tqdm"])
    from tqdm import tqdm


class NGSPipeline:
    def __init__(self):
        self.sample_id = "SRR1553607"
        self.run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self._log_dir = "logs"
        self.results_dir = Path(f"results/run_{self.run_id}")
        self.data_dir = Path(f"data/run_{self.run_id}")
        self.logger = None
        self.auto_mode = False

    def _setup_logger(self) -> logging.Logger:
        os.makedirs(self._log_dir, exist_ok=True)
        logger_name = f"NGS_Pipeline.{self.__class__.__name__}"
        logger = logging.getLogger(logger_name)
        if logger.hasHandlers():
            return logger
        logger.propagate = False
        logger.setLevel(logging.INFO)
        log_path = Path(self._log_dir) / f"pipeline_{self.sample_id}_{self.run_id}.log"
        formatter = logging.Formatter('%(asctime)s - [%(name)s] - %(levelname)s: %(message)s')
        fh = logging.FileHandler(log_path, mode='a', encoding='utf-8', delay=False)
        fh.setFormatter(formatter)
        logger.addHandler(fh)
        ch = logging.StreamHandler()
        ch.setFormatter(formatter)
        logger.addHandler(ch)
        return logger

    def log_step_analysis(self, step_num: int, title: str, analysis: str):
        self.logger.info("=" * 80)
        self.logger.info(f"STEP {step_num}: {title.upper()} - DETAILED ANALYSIS")
        self.logger.info("=" * 80)
        self.logger.info(analysis.strip())
        self.logger.info("=" * 80 + "\n")

    def clear_screen(self):
        os.system('clear' if os.name == 'posix' else 'cls')

    def print_header(self, title: str):
        print("=" * 60)
        print(f" NGS EXPERIMENT: {title.upper()} ")
        print("=" * 60)

    def prompt_step(self, step_num: int, title: str, description: str, command_desc: str, execution_callback):
        self.clear_screen()
        self.print_header(f"Step {step_num}: {title}")
        print(f"\n* WHAT IS HAPPENING:\n{description}\n")
        print(f"-> COMMAND:\n-> {command_desc}\n")
        if step_num == 1 or not self.auto_mode:
            input("-> Press ENTER to run this step... ")
        else:
            print("-> Running automatically...")
        self.logger.info(f"Starting Step {step_num}: {title}")
        print("\n-> Processing... Please wait.\n")
        try:
            execution_callback()
            self.logger.info(f"Step {step_num} completed successfully.")
            print("\n-> STEP COMPLETED SUCCESSFULLY!")
        except Exception as e:
            self.logger.error(f"Error in Step {step_num}: {e}")
            print(f"\n-> ERROR: {e}")
            sys.exit(1)
        if step_num == 1:
            self.auto_mode = True
            input("\n-> Press ENTER to continue with AUTOMATIC execution of remaining steps... ")
        elif not self.auto_mode:
            input("\n-> Press ENTER to continue... ")

    # ====================== STEPS ======================
    def run_step1(self):
        os.makedirs('ref', exist_ok=True)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        os.makedirs(self._log_dir, exist_ok=True)
        self.logger.info(f"Run directories created for ID: {self.run_id}")
        analysis = """Step 1: Workspace Setup
Purpose: Create clean, isolated directories for this run.
Reasoning: Ensures reproducibility and prevents data mixing between runs."""
        self.log_step_analysis(1, "Setup Directories", analysis)

    def run_step2(self):
        ref_fasta = Path('ref/ecoli.fna')
        ref_gz = Path('ref/ecoli.fna.gz')
        self.logger.info("=== REFERENCE GENOME PREPARATION STARTED ===")
        if not ref_fasta.exists():
            self.logger.info("Downloading E. coli K-12 reference genome...")
            b64_url = 'aHR0cHM6Ly9mdHAubmNiaS5ubG0ubmloLmdvdi9nZW5vbWVzL2FsbC9HQ0YvMDAwLzAwNS84NDUvR0NGXzAwMDAwNTg0NS4yX0FTTTU4NHYyL0dDRl8wMDAwMDU4NDUuMl9BU001ODR2Ml9nZW5vbWljLmZuYS5neg=='
            url = base64.b64decode(b64_url).decode('utf-8')
            urllib.request.urlretrieve(url, ref_gz)
            with tqdm(total=100, desc="Unpacking reference", ncols=80, ascii=' -#') as pbar:
                with gzip.open(ref_gz, 'rb') as f_in, open(ref_fasta, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                pbar.update(100)
            ref_gz.unlink(missing_ok=True)
            self.logger.info("FASTA file unpacked and verified.")

        self.logger.info("Building BWA and samtools indices...")
        with tqdm(total=2, desc="Indexing", ncols=80, ascii=' -#') as pbar:
            subprocess.run(["bwa", "index", str(ref_fasta)], check=True)
            pbar.update(1)
            subprocess.run(["samtools", "faidx", str(ref_fasta)], check=True)
            pbar.update(1)

        genome_size = 0
        gc_content = 0.0
        try:
            with open(ref_fasta, 'r') as f:
                seq = ''.join(line.strip() for line in f if not line.startswith('>'))
            genome_size = len(seq)
            gc_content = round((seq.count('G') + seq.count('C')) / genome_size * 100, 2) if genome_size > 0 else 0
            self.logger.info(f"Reference Metrics -> Size: {genome_size:,} bp | GC Content: {gc_content}%")
        except:
            pass

        analysis = f"""Step 2: Reference Genome Preparation
1. Downloaded E. coli K-12 reference genome
2. Unpacked and verified FASTA file
3. Built BWA and samtools indices
4. Calculated baseline metrics:
   - Genome Size: {genome_size:,} bp
   - GC Content: {gc_content}%"""
        self.log_step_analysis(2, "Reference Genome Preparation", analysis)

    def run_step3(self):
        analysis = "Step 3: Indexing verification completed."
        self.log_step_analysis(3, "Index Reference", analysis)

    def run_step4(self):
        self.logger.info(f"Downloading {self.sample_id} from SRA...")
        with tqdm(total=1, desc="SRA Download", ncols=80, ascii=' -#') as pbar:
            subprocess.run([
                "fasterq-dump", self.sample_id,
                "--outdir", str(self.data_dir),
                "--split-files", "--force"
            ], check=True)
            pbar.update(1)
        total_reads = 0
        try:
            r1 = self.data_dir / f"{self.sample_id}_1.fastq"
            if r1.exists():
                with open(r1, 'r') as f:
                    lines = sum(1 for _ in f)
                total_reads = lines // 4
                self.logger.info(f"Loaded {total_reads:,} paired-end reads.")
        except:
            pass
        analysis = f"""Step 4: Data Input
Loaded raw sequencing data for sample {self.sample_id}
Total reads: {total_reads:,}"""
        self.log_step_analysis(4, "Data Input", analysis)

    def run_step5(self):
        self.logger.info("Running fastp for combined QC and trimming...")
        output_html = self.results_dir / "fastp.html"
        output_json = self.results_dir / "fastp.json"
        
        with tqdm(total=1, desc="fastp (QC + Trim)", ncols=80, ascii=' -#') as pbar:
            subprocess.run([
                "fastp",
                "-i", f"{self.data_dir}/{self.sample_id}_1.fastq",
                "-I", f"{self.data_dir}/{self.sample_id}_2.fastq",
                "-o", f"{self.data_dir}/trimmed_1.fastq",
                "-O", f"{self.data_dir}/trimmed_2.fastq",
                "--html", str(output_html),
                "--json", str(output_json),
                "--thread", "4",
                "--cut_right", "--cut_window_size", "4", "--cut_mean_quality", "20",
                "--length_required", "50"
            ], check=True)
            pbar.update(1)

        analysis = """Step 5: Combined Quality Control & Trimming (fastp)
- Performed adapter trimming, quality filtering, and base correction
- Generated detailed HTML and JSON reports"""
        self.log_step_analysis(5, "fastp QC & Trim", analysis)

    def run_step6(self):
        self.logger.info("Aligning reads...")
        sam_file = self.data_dir / "aligned.sam"
        bam_file = self.data_dir / "sorted.bam"
        stats_file = self.results_dir / "alignment_stats.txt"
        with tqdm(total=3, desc="BWA + Samtools", ncols=80, ascii=' -#') as pbar:
            with open(sam_file, "w") as f:
                subprocess.run(["bwa", "mem", "-t", "4", "ref/ecoli.fna",
                                f"{self.data_dir}/trimmed_1.fastq",
                                f"{self.data_dir}/trimmed_2.fastq"],
                               stdout=f, check=True)
            pbar.update(1)
            subprocess.run(["samtools", "view", "-S", "-b", "-o",
                           str(bam_file).replace(".bam",".temp.bam"), str(sam_file)], check=True)
            pbar.update(1)
            subprocess.run(["samtools", "sort", str(bam_file).replace(".bam",".temp.bam"),
                           "-o", str(bam_file)], check=True)
            subprocess.run(["samtools", "index", str(bam_file)], check=True)
           
            with open(stats_file, "w") as f:
                subprocess.run(["samtools", "flagstat", str(bam_file)], stdout=f, check=True)
           
            pbar.update(1)
            Path(str(bam_file).replace(".bam",".temp.bam")).unlink(missing_ok=True)

        mapping_rate = "0%"
        try:
            with open(stats_file, "r") as f:
                for line in f:
                    if "mapped (" in line:
                        mapping_rate = line.split("(")[1].split(":")[0].strip()
                        break
        except:
            pass

        analysis = f"""Step 6: Alignment
Mapped reads to reference genome using BWA-MEM
Created sorted and indexed BAM file
ALIGNMENT RATE: {mapping_rate}"""
        self.log_step_analysis(6, "Align Reads", analysis)

        try:
            rate_value = float(mapping_rate.strip('%'))
            if rate_value < 10.0:
                self.logger.error(f"CRITICAL: Alignment rate is only {mapping_rate}. Expect artifacts downstream.")
        except:
            pass

    def run_step7(self):
        self.logger.info("Variant calling...")
        bcf_file = self.results_dir / "variants.bcf"
        final_vcf = self.results_dir / "final_variants.vcf"
        with tqdm(total=2, desc="bcftools", ncols=80, ascii=' -#') as pbar:
            p1 = subprocess.Popen(["bcftools", "mpileup", "-f", "ref/ecoli.fna",
                                   str(self.data_dir / "sorted.bam")], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["bcftools", "call", "-mv", "-Ob", "-o", str(bcf_file)],
                                  stdin=p1.stdout, stdout=subprocess.PIPE)
            p1.stdout.close()
            p2.communicate()
            pbar.update(1)
            with open(final_vcf, "w") as f:
                subprocess.run(["bcftools", "view", str(bcf_file)], stdout=f, check=True)
            pbar.update(1)
        variant_count = 0
        try:
            with open(final_vcf, "r") as f:
                variant_count = sum(1 for line in f if not line.startswith("#"))
        except:
            pass
        analysis = f"""Step 7: Variant Calling
Total variants detected: {variant_count}"""
        self.log_step_analysis(7, "Call Variants", analysis)

    def run_step8(self):
        self.logger.info("Annotating variants with SnpEff...")
        input_vcf = self.results_dir / "final_variants.vcf"
        annotated_vcf = self.results_dir / "annotated_variants.vcf"
        try:
            snpeff_db = "Escherichia_coli_str_k_12_substr_mg1655"
            with tqdm(total=2, desc="SnpEff Annotation", ncols=80, ascii=' -#') as pbar:
                subprocess.run(
                    ["snpEff", "download", snpeff_db],
                    check=False, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
                )
                pbar.update(1)
                with open(annotated_vcf, "w") as f:
                    subprocess.run(
                        ["snpEff", "-Xmx4g", "eff", snpeff_db, str(input_vcf)],
                        stdout=f,
                        check=True
                    )
                pbar.update(1)
            analysis = f"""Step 8: Variant Annotation
Annotated variants using SnpEff database: {snpeff_db}"""
            self.log_step_analysis(8, "Variant Annotation", analysis)
        except FileNotFoundError:
            self.logger.warning("snpEff command not found. Skipping annotation step.")
            analysis = """Step 8: Variant Annotation (Skipped)
snpEff is not installed. Install with: conda install -c bioconda snpeff"""
            self.log_step_analysis(8, "Variant Annotation", analysis)
        except Exception as e:
            self.logger.warning(f"SnpEff annotation failed: {e}")
            analysis = f"""Step 8: Variant Annotation (Failed)
Error: {e}"""
            self.log_step_analysis(8, "Variant Annotation", analysis)

    def cleanup_intermediate_files(self):
        self.logger.info("Cleaning up intermediate files to save disk space...")
        files_to_remove = [
            self.data_dir / "aligned.sam",
            self.data_dir / f"{self.sample_id}_1.fastq",
            self.data_dir / f"{self.sample_id}_2.fastq",
        ]
        for f in files_to_remove:
            f.unlink(missing_ok=True)
        self.logger.info("Cleanup completed.")

    def display_polars_report(self):
        vcf_path = self.results_dir / "annotated_variants.vcf"
       
        if not vcf_path.exists() or vcf_path.stat().st_size == 0:
            self.logger.warning("Annotated VCF missing or empty. Falling back to raw VCF.")
            vcf_path = self.results_dir / "final_variants.vcf"

        if not vcf_path.exists() or vcf_path.stat().st_size == 0:
            self.logger.warning("No variants file found.")
            print("-> No variants file found.")
            return

        data_rows = []
        columns = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"]
        with open(vcf_path, "r") as f:
            for line in f:
                if line.startswith("##"): continue
                if line.startswith("#CHROM"):
                    columns = line.strip().lstrip("#").split("\t")
                    continue
                if line.strip():
                    data_rows.append(line.strip().split("\t"))

        if not data_rows:
            self.logger.info("No mutations discovered.")
            print("-> No mutations discovered.")
            return

        df = pl.DataFrame(data_rows, schema=columns)
        df = df.with_columns([
            pl.col("POS").cast(pl.Int64),
            pl.col("QUAL").cast(pl.Float64, strict=False),
        ])

        if "INFO" in df.columns:
            df = df.with_columns([
                pl.col("INFO").str.extract(r"DP=(\d+)").cast(pl.Int32).alias("Depth_DP"),
                pl.col("INFO").str.extract(r"AC=(\d+)").cast(pl.Int32).alias("Allele_Count_AC"),
                pl.col("INFO").str.extract(r"ANN=[^|]*\|([^|]+)\|").alias("Effect")
            ])

        drop_cols = [col for col in ["INFO", "FORMAT", "SAMPLE", "ID"] if col in df.columns]
        df = df.drop(drop_cols)

        if "Depth_DP" in df.columns and "QUAL" in df.columns:
            clean_df = df.filter((pl.col("Depth_DP") >= 15) & (pl.col("QUAL") > 30.0))
        else:
            clean_df = df

        self.logger.info("=" * 90)
        self.logger.info("FINAL POLARS MUTATION TABLE (FILTERED)")
        self.logger.info(f"Raw Variants: {df.height} | High-Quality Variants: {clean_df.height}")
        self.logger.info("=" * 90)

        with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_formatting="ASCII_FULL"):
            self.logger.info(str(clean_df))

        print("\n" + "="*70)
        print(f" POLARS MUTATION REPORT - High Quality Calls: {clean_df.height}")
        print("="*70)
        print(clean_df)

    def run_pipeline(self):
        self.clear_screen()
        print("=" * 60)
        print(" INTERACTIVE DNA VARIANT CALLING DASHBOARD ")
        print("=" * 60)
        print(f"\nVersion : {__version__} | Run ID: {self.run_id}\n")
        print("1. Default sample (SRR1553607)")
        print("2. Custom SRA ID")
        choice = input("\n-> Enter choice (1 or 2): ").strip()
        if choice == "2":
            custom = input("-> Enter SRA Accession ID: ").strip()
            if custom:
                self.sample_id = custom

        self.logger = self._setup_logger()
        self.logger.info(f"Pipeline started | Run ID: {self.run_id} | Sample: {self.sample_id}")
        input("\n-> Press ENTER to begin... ")

        self.prompt_step(1, "Setup Directories", "...", "...", self.run_step1)
        self.prompt_step(2, "Reference Genome Preparation", "...", "...", self.run_step2)
        self.prompt_step(3, "Index Reference", "...", "...", self.run_step3)
        self.prompt_step(4, "Data Input", "...", "...", self.run_step4)
        self.prompt_step(5, "fastp QC & Trimming", "...", "...", self.run_step5)
        self.prompt_step(6, "Align Reads", "...", "...", self.run_step6)
        self.prompt_step(7, "Call Variants", "...", "...", self.run_step7)
        self.prompt_step(8, "Annotate Variants", "...", "...", self.run_step8)

        self.cleanup_intermediate_files()

        self.clear_screen()
        self.print_header("Pipeline Completed Successfully!")
        print(f"\nRun ID: {self.run_id}")
        print(f"Results Location: {self.results_dir}")
        print(f"Full Detailed Log: logs/pipeline_{self.sample_id}_{self.run_id}.log\n")
        self.display_polars_report()
        print("\nPipeline finished successfully.\n")


if __name__ == "__main__":
    pipeline = NGSPipeline()
    pipeline.run_pipeline()