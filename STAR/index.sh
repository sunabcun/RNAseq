#!/bin/bash
#SBATCH --job-name="trap"
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH -p medicine_q
#SBATCH --mail-type="ALL"
#SBATCH -t 48:48:15
#SBATCH -J starindex

STAR --runThreadN 5 \
--runMode genomeGenerate \
--genomeDir ./ \
--genomeFastaFiles \
GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile gencode.v34.primary_assembly.annotation.gtf \
--sjdbOverhang 99
