#!/bin/bash
set -e

echo "==> Creating project directories..."
mkdir -p data ref results

if [ ! -f ref/ecoli.fna ]; then
    echo "==> Fetching E. coli K-12 MG1655 genome..."
    rm -rf ref && mkdir -p ref
    
    wget -O ref/ecoli.fna "https://ebi.ac.uk"
    
    echo "==> Indexing reference genome..."
    bwa index ref/ecoli.fna
    samtools faidx ref/ecoli.fna
fi

echo "==> Downloading SRA dataset SRR1553607..."
fasterq-dump SRR1553607 --outdir data --split-files

echo "==> Running FastQC..."
mkdir -p results/fastqc
fastqc data/SRR1553607_1.fastq data/SRR1553607_2.fastq -o results/fastqc/

echo "==> Quality trimming with Trimmomatic..."
trimmomatic PE -phred33 \
  data/SRR1553607_1.fastq data/SRR1553607_2.fastq \
  data/trimmed_1.fastq data/unpaired_1.fastq \
  data/trimmed_2.fastq data/unpaired_2.fastq \
  SLIDINGWINDOW:4:20 MINLEN:50

echo "==> Aligning reads to reference with BWA-MEM..."
bwa mem -t 4 ref/ecoli.fna data/trimmed_1.fastq data/trimmed_2.fastq > data/aligned.sam

samtools view -S -b data/aligned.sam > data/aligned.bam
samtools sort data/aligned.bam -o data/sorted.bam
samtools index data/sorted.bam

echo "==> Calling variants with bcftools..."
bcftools mpileup -f ref/ecoli.fna data/sorted.bam | bcftools call -mv -Ob -o results/variants.bcf
bcftools view results/variants.bcf > results/final_variants.vcf

echo "==> Pipeline complete! Final variants are in results/final_variants.vcf"
