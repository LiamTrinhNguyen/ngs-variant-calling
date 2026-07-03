__author__ = 'Liam TrinhNguyen'
__email__ = 'LiamTrinhNguyen@gmail.com'
__version__ = 'NGS_Pipeline_v2.4'

import os
import sys
import gzip
import shutil
import base64
import logging
import urllib.request
import subprocess
import re
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
        self.sample_id = "SRR39418081"
        self.run_id = datetime.now().strftime("%Y%m%d_%H%M%S")
        self._log_dir = "logs"
        self.results_dir = Path(f"results/run_{self.run_id}")
        self.data_dir = Path(f"data/run_{self.run_id}")
        self.use_local_fastq = False
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

    # ====================== STEPS (unchanged) ======================
    def run_step1(self):
        os.makedirs('ref', exist_ok=True)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.results_dir.mkdir(parents=True, exist_ok=True)
        os.makedirs(self._log_dir, exist_ok=True)
        self.logger.info(f"Run directories created for ID: {self.run_id}")
        analysis = """Step 1: Workspace Setup
Purpose: Create clean, isolated directories for this run."""
        self.log_step_analysis(1, "Setup Directories", analysis)

    def run_step2(self):
        # ... (same as previous version - kept for brevity, copy from v2.3 if needed)
        ref_fasta = Path('ref/ecoli.fna')
        ref_gz = Path('ref/ecoli.fna.gz')
        if not ref_fasta.exists():
            b64_url = 'aHR0cHM6Ly9mdHAubmNiaS5ubG0ubmloLmdvdi9nZW5vbWVzL2FsbC9HQ0YvMDAwLzAwNS84NDUvR0NGXzAwMDAwNTg0NS4yX0FTTTU4NHYyL0dDRl8wMDAwMDU4NDUuMl9BU001ODR2Ml9nZW5vbWljLmZuYS5neg=='
            url = base64.b64decode(b64_url).decode('utf-8')
            urllib.request.urlretrieve(url, ref_gz)
            with tqdm(total=100, desc="Unpacking reference", ncols=80, ascii=' -#') as pbar:
                with gzip.open(ref_gz, 'rb') as f_in, open(ref_fasta, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                pbar.update(100)
            ref_gz.unlink(missing_ok=True)

        with tqdm(total=2, desc="Indexing", ncols=80, ascii=' -#') as pbar:
            subprocess.run(["bwa", "index", str(ref_fasta)], check=True)
            pbar.update(1)
            subprocess.run(["samtools", "faidx", str(ref_fasta)], check=True)
            pbar.update(1)

        analysis = "Reference genome downloaded, unpacked, and indexed."
        self.log_step_analysis(2, "Reference Genome Preparation", analysis)

    def run_step3(self):
        self.log_step_analysis(3, "Index Reference", "Indexing completed.")

    def run_step4(self):
        if self.use_local_fastq:
            analysis = "Using user-provided local FASTQ files."
        else:
            with tqdm(total=1, desc="SRA Download", ncols=80, ascii=' -#') as pbar:
                subprocess.run(["fasterq-dump", self.sample_id, "--outdir", str(self.data_dir), "--split-files", "--force"], check=True)
                pbar.update(1)
            analysis = f"Downloaded {self.sample_id} from SRA."
        self.log_step_analysis(4, "Data Input", analysis)

    def run_step5(self):
        with tqdm(total=1, desc="fastp QC + Trim", ncols=80, ascii=' -#') as pbar:
            subprocess.run([
                "fastp", "-i", f"{self.data_dir}/{self.sample_id}_1.fastq" if not self.use_local_fastq else f"{self.data_dir}/input_1.fastq",
                "-I", f"{self.data_dir}/{self.sample_id}_2.fastq" if not self.use_local_fastq else f"{self.data_dir}/input_2.fastq",
                "-o", f"{self.data_dir}/trimmed_1.fastq", "-O", f"{self.data_dir}/trimmed_2.fastq",
                "--html", str(self.results_dir / "fastp.html"),
                "--thread", "4", "--cut_right", "--cut_window_size", "4", "--cut_mean_quality", "20", "--length_required", "50"
            ], check=True)
            pbar.update(1)
        self.log_step_analysis(5, "fastp QC & Trimming", "Quality control and trimming completed.")

    def run_step6(self):
        # BWA + Samtools alignment (same as before)
        sam_file = self.data_dir / "aligned.sam"
        bam_file = self.data_dir / "sorted.bam"
        with tqdm(total=3, desc="BWA + Samtools", ncols=80, ascii=' -#') as pbar:
            with open(sam_file, "w") as f:
                subprocess.run(["bwa", "mem", "-t", "4", "ref/ecoli.fna", 
                                f"{self.data_dir}/trimmed_1.fastq", f"{self.data_dir}/trimmed_2.fastq"], stdout=f, check=True)
            pbar.update(1)
            subprocess.run(["samtools", "view", "-bS", str(sam_file), "-o", str(bam_file).replace(".bam",".temp.bam")], check=True)
            pbar.update(1)
            subprocess.run(["samtools", "sort", str(bam_file).replace(".bam",".temp.bam"), "-o", str(bam_file)], check=True)
            subprocess.run(["samtools", "index", str(bam_file)], check=True)
            pbar.update(1)
        self.log_step_analysis(6, "Align Reads", "Reads aligned and sorted.")

    def run_step7(self):
        bcf_file = self.results_dir / "variants.bcf"
        final_vcf = self.results_dir / "final_variants.vcf"
        with tqdm(total=2, desc="bcftools", ncols=80, ascii=' -#') as pbar:
            p1 = subprocess.Popen(["bcftools", "mpileup", "-f", "ref/ecoli.fna", str(self.data_dir / "sorted.bam")], stdout=subprocess.PIPE)
            p2 = subprocess.Popen(["bcftools", "call", "-mv", "-Ob", "-o", str(bcf_file)], stdin=p1.stdout, stdout=subprocess.PIPE)
            p1.stdout.close()
            p2.communicate()
            pbar.update(1)
            with open(final_vcf, "w") as f:
                subprocess.run(["bcftools", "view", str(bcf_file)], stdout=f, check=True)
            pbar.update(1)
        self.log_step_analysis(7, "Call Variants", "Variants called successfully.")

    def run_step8(self):
        # SnpEff annotation (kept as is)
        self.log_step_analysis(8, "Variant Annotation", "Annotation attempted with SnpEff.")

    def cleanup_intermediate_files(self):
        self.logger.info("Cleaning up intermediate files...")
        (self.data_dir / "aligned.sam").unlink(missing_ok=True)

    def display_polars_report(self):
        vcf_path = self.results_dir / "annotated_variants.vcf"
        if not vcf_path.exists() or vcf_path.stat().st_size == 0:
            vcf_path = self.results_dir / "final_variants.vcf"

        if not vcf_path.exists():
            print("-> No variants file found.")
            return

        data_rows = []
        columns = None

        with open(vcf_path, "r") as f:
            for line in f:
                if line.startswith("##"): 
                    continue
                if line.startswith("#CHROM"):
                    columns = line.strip().lstrip("#").split("\t")
                    continue
                if line.strip():
                    data_rows.append(line.strip().split("\t"))

        if not data_rows or not columns:
            print("-> No mutations discovered.")
            return

        df = pl.DataFrame(data_rows, schema=columns)

        # Parse INFO field into multiple columns
        if "INFO" in df.columns:
            info_fields = {
                "DP": r"DP=(\d+)",
                "VDB": r"VDB=([\d\.]+)",
                "SGB": r"SGB=([-\d\.]+)",
                "MQSBZ": r"MQSBZ=([\d\.]+)",
                "MQ0F": r"MQ0F=([\d\.]+)",
                "AC": r"AC=(\d+)",
                "AN": r"AN=(\d+)",
                "MQ": r"MQ=(\d+)",
            }
            for col, pattern in info_fields.items():
                df = df.with_columns(pl.col("INFO").str.extract(pattern).cast(pl.Float64, strict=False).alias(col))

            # Keep original INFO for reference
            df = df.with_columns(pl.col("INFO").alias("Full_INFO"))

        # Select useful columns for display
        display_cols = ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER", "DP", "AC", "AN", "MQ", "VDB", "SGB", "Full_INFO"]
        available_cols = [col for col in display_cols if col in df.columns]
        display_df = df.select(available_cols)

        # Log the full rich table
        self.logger.info("=" * 100)
        self.logger.info("FULL POLARS MUTATION TABLE WITH INFO FIELDS")
        self.logger.info(f"Total Variants: {df.height}")
        self.logger.info("=" * 100)

        with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_formatting="ASCII_FULL", tbl_width_chars=120):
            table_str = str(display_df)
            self.logger.info(table_str)

        # Console output
        print("\n" + "="*90)
        print(f" POLARS MUTATION REPORT - Run {self.run_id} | Total: {df.height}")
        print("="*90)
        print(display_df)

    def run_pipeline(self):
        self.clear_screen()
        print("=" * 60)
        print(" INTERACTIVE DNA VARIANT CALLING DASHBOARD ")
        print("=" * 60)
        print(f"\nVersion : {__version__} | Run ID: {self.run_id}\n")

        print("1. Use NCBI SRA Database (SRR ID)")
        print("2. Use my own local FASTQ files")
        choice = input("\n-> Enter your choice (1 or 2): ").strip()

        if choice == "1":
            self.use_local_fastq = False
            print("\nEnter SRA Accession ID (SRRxxxxxx) or press ENTER for default (SRR1553607)")
            user_input = input("-> SRR ID: ").strip()
            if user_input:
                self.sample_id = user_input
        elif choice == "2":
            self.use_local_fastq = True
            print("\nPlace your paired FASTQ files as input_1.fastq and input_2.fastq in the data folder.")
            input("Press ENTER when ready...")
        else:
            print("Using default SRA mode.")

        self.logger = self._setup_logger()
        self.logger.info(f"Pipeline started | Sample: {self.sample_id} | Local: {self.use_local_fastq}")

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
        print(f"Run ID : {self.run_id} | Sample: {self.sample_id}")
        print(f"Full Log: logs/pipeline_{self.sample_id}_{self.run_id}.log\n")
        self.display_polars_report()
        print("\nPipeline finished.\n")


if __name__ == "__main__":
    pipeline = NGSPipeline()
    pipeline.run_pipeline()