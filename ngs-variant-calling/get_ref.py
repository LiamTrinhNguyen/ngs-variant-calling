import urllib.request
import gzip
import shutil
import os

# Base paths broken into chunks to avoid any local clipboard filtering bugs
url_host = "https://ftp.ncbi.nlm.nih.gov"
url_path = "/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/"
url_file = "GCF_000005845.2_ASM584v2_genomic.fna.gz"
full_url = url_host + url_path + url_file

os.makedirs("ref", exist_ok=True)
gz_target = "ref/ecoli.fna.gz"
fna_target = "ref/ecoli.fna"

print("==> Downloading reference genome from NCBI RefSeq...")
urllib.request.urlretrieve(full_url, gz_target)

print("==> Uncompressing reference genome...")
with gzip.open(gz_target, 'rb') as f_in:
    with open(fna_target, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

# Clean up compressed file archive
if os.path.exists(gz_target):
    os.remove(gz_target)

print("==> Success! Reference genome saved to ref/ecoli.fna")
# 1. Index the downloaded sequence structure
bwa index ref/ecoli.fna
samtools faidx ref/ecoli.fna

# 2. Map trimmed sequence reads to the reference file
bwa mem -t 4 ref/ecoli.fna data/trimmed_1.fastq data/trimmed_2.fastq > data/aligned.sam

# 3. Process maps into indexed BAM format
samtools view -S -b data/aligned.sam > data/aligned.bam
samtools sort data/aligned.bam -o data/sorted.bam
samtools index data/sorted.bam

# 4. Execute final variant calling calculation
bcftools mpileup -f ref/ecoli.fna data/sorted.bam | bcftools call -mv -Ob -o results/variants.bcf
bcftools view results/variants.bcf > results/final_variants.vcf

echo "==> Pipeline complete! Final variants are in results/final_variants.vcf"
