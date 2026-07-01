#!/bin/bash
set -e
mkdir -p ref

echo "==> Pulling E. coli reference genome..."
wget -O ref/ecoli.fna "https://ebi.ac.uk"

echo "==> Building BWA Index..."
bwa index ref/ecoli.fna

echo "==> Building Samtools Index..."
samtools faidx ref/ecoli.fna

echo "==> Reference Ready!"
