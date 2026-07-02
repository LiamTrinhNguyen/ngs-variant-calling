__author__ = 'Liam TrinhNguyen'
__email__ = 'LiamTrinhNguyen@gmail.com'
__version__ = 'NGS_Pipeline_v1.6'

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
        """Save detailed analysis to log file"""
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
Purpose: Create clean, isolated directories for this specific run.
Reasoning: Ensures full reproducibility and prevents data contamination between runs.
Best Practice: Standard in professional bioinformatics workflows."""
        self.log_step_analysis(1, "Setup Directories", analysis)

    def run_step2(self):
        ref_file = Path('ref/ecoli.fna')
        if ref_file.exists():
            self.logger.info("Reference genome already exists.")
        else:
            self.logger.info("Downloading reference genome...")
            b64_url = 'aHR0cHM6Ly9mdHAubmNiaS5ubG0ubmloLmdvdi9nZW5vbWVzL2FsbC9HQ0YvMDAwLzAwNS84NDUvR0NGXzAwMDAwNTg0NS4yX0FTTTU4NHYyL0dDRl8wMDAwMDU4NDUuMl9BU001ODR2Ml9nZW5vbWljLmZuYS5neg=='
            url = base64.b64decode(b64_url).decode('utf-8')
            urllib.request.urlretrieve(url, 'ref/ecoli.fna.gz')
            
            with tqdm(total=100, desc="Unpacking reference", ncols=80, ascii=' -#') as pbar:
                with gzip.open('ref/ecoli.fna.gz', 'rb') as f_in, open(ref_file, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
                pbar.update(100)
            
            Path('ref/ecoli.fna.gz').unlink(missing_ok=True)
            self.logger.info("Reference genome ready.")

        analysis = """Step 2: Reference Genome Download
Purpose: Obtain high-quality E. coli K-12 reference genome.
Biological Context: Serves as the trusted standard for variant comparison.
Reasoning: All detected mutations are relative to this reference strain."""
        self.log_step_analysis(2, "Download Reference Genome", analysis)

    def run_step3(self):
        self.logger.info("Indexing reference genome...")
        with tqdm(total=2, desc="Indexing", ncols=80, ascii=' -#') as pbar:
            subprocess.run(["bwa", "index", "ref/ecoli.fna"], check=True)
            pbar.update(1)
            subprocess.run(["samtools", "faidx", "ref/ecoli.fna"], check=True)
            pbar.update(1)

        analysis = """Step 3: Reference Indexing
Purpose: Build searchable index files for fast alignment.
Reasoning: Dramatically improves speed and efficiency of read mapping."""
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

        analysis = """Step 4: Sequencing Data Input
Purpose: Acquire raw FASTQ reads for analysis.
Context: Public SRA dataset or user-provided local FASTQ files.
Reasoning: This is the experimental data containing potential biological variants."""
        self.log_step_analysis(4, "Download Reads", analysis)

    def run_step5(self):
        qc_dir = self.results_dir / "fastqc"
        qc_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info("Running FastQC...")
        with tqdm(total=1, desc="FastQC", ncols=80, ascii=' -#') as pbar:
            subprocess.run([
                "fastqc",
                f"{self.data_dir}/{self.sample_id}_1.fastq",
                f"{self.data_dir}/{self.sample_id}_2.fastq",
                "-o", str(qc_dir)
            ], check=True)
            pbar.update(1)

        analysis = """Step 5: Quality Control (FastQC)
Purpose: Assess quality of raw sequencing reads.
Key Checks: Per-base quality, GC content, adapters, duplication levels.
Reasoning: Early detection of problems prevents errors in downstream analysis."""
        self.log_step_analysis(5, "FastQC Quality Check", analysis)

    def run_step6(self):
        self.logger.info("Trimming reads...")
        with tqdm(total=1, desc="Trimmomatic", ncols=80, ascii=' -#') as pbar:
            subprocess.run([
                "trimmomatic", "PE", "-phred33",
                f"{self.data_dir}/{self.sample_id}_1.fastq",
                f"{self.data_dir}/{self.sample_id}_2.fastq",
                f"{self.data_dir}/trimmed_1.fastq", f"{self.data_dir}/unpaired_1.fastq",
                f"{self.data_dir}/trimmed_2.fastq", f"{self.data_dir}/unpaired_2.fastq",
                "SLIDINGWINDOW:4:20", "MINLEN:50"
            ], check=True)
            pbar.update(1)

        analysis = """Step 6: Read Trimming (Trimmomatic)
Purpose: Remove low-quality bases and adapter sequences.
Reasoning: Improves alignment accuracy and reduces false positive variant calls."""
        self.log_step_analysis(6, "Trim Reads", analysis)

    def run_step7(self):
        self.logger.info("Aligning reads...")
        sam_file = self.data_dir / "aligned.sam"
        bam_file = self.data_dir / "sorted.bam"

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
            pbar.update(1)

            Path(str(bam_file).replace(".bam",".temp.bam")).unlink(missing_ok=True)

        analysis = """Step 7: Alignment (BWA-MEM + Samtools)
Purpose: Map cleaned reads to the reference genome.
Reasoning: Accurate mapping is the foundation for reliable variant detection."""
        self.log_step_analysis(7, "Align Reads", analysis)

    def run_step8(self):
        self.logger.info("Variant calling...")
        bcf_file = self.results_dir / "variants.bcf"
        final_vcf = self.results_dir / "final_variants.vcf"

        with tqdm(total=2, desc="bcftools mpileup + call", ncols=80, ascii=' -#') as pbar:
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

        analysis = """Step 8: Variant Calling (BCFtools)
Purpose: Identify genomic differences between sample and reference.
Reasoning: Statistical approach distinguishes true mutations from sequencing errors."""
        self.log_step_analysis(8, "Call Variants", analysis)

    def display_polars_report(self):
        vcf_path = self.results_dir / "final_variants.vcf"
        if not vcf_path.exists():
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
        if "POS" in df.columns:
            df = df.with_columns(pl.col("POS").cast(pl.Int64))

        self.logger.info("=" * 90)
        self.logger.info("FINAL POLARS MUTATION TABLE")
        self.logger.info(f"Total Variants: {df.height}")
        self.logger.info("=" * 90)

        with pl.Config(tbl_rows=-1, tbl_cols=-1, tbl_formatting="ASCII_FULL"):
            self.logger.info(str(df))

        self.logger.info("=" * 90)

        print("\n" + "="*70)
        print(f" POLARS MUTATION REPORT - Run {self.run_id} (Total Variants: {df.height})")
        print("="*70)

        display_cols = [c for c in ["CHROM", "POS", "REF", "ALT", "QUAL", "FILTER"] if c in df.columns]
        with pl.Config(tbl_rows=-1, tbl_cols=-1):
            print(df.select(display_cols))

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
        self.prompt_step(2, "Download Reference Genome", "...", "...", self.run_step2)
        self.prompt_step(3, "Index Reference", "...", "...", self.run_step3)
        self.prompt_step(4, "Download Reads", "...", "...", self.run_step4)
        self.prompt_step(5, "FastQC Quality Check", "...", "...", self.run_step5)
        self.prompt_step(6, "Trim Reads", "...", "...", self.run_step6)
        self.prompt_step(7, "Align Reads", "...", "...", self.run_step7)
        self.prompt_step(8, "Call Variants", "...", "...", self.run_step8)

        self.clear_screen()
        self.print_header("Pipeline Completed Successfully!")
        print(f"\nRun ID: {self.run_id}")
        print(f"Results Location: {self.results_dir}")
        print(f"Full Detailed Log (with analysis): logs/pipeline_{self.sample_id}_{self.run_id}.log\n")

        self.display_polars_report()
        print("\n🎉 Pipeline finished successfully.\n")


if __name__ == "__main__":
    pipeline = NGSPipeline()
    pipeline.run_pipeline()