#!/bin/bash
#SBATCH --job-name="trap"
#SBATCH -N 6
#SBATCH --mem=16G
#SBATCH -p medicine_q
#SBATCH --mail-type="ALL"
#SBATCH -t 48:48:15
#SBATCH -J HISAT2align

for i in *_R1.fastq.gz; do
/gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Tools/hisat2-2.2.1/hisat2 \
-p 8 \
--rg-id=${i%_R1.fastq.gz} \
--rg SM:${i%_R1.fastq.gz}_SM \
--rg LB:${i%_R1.fastq.gz}_LB \
--rg PL:ILLUMINA \
-x /gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Reference/grch38/genome \
--dta \
--rna-strandness RF \
-1 $i \
-2 ${i%_R1.fastq.gz}_R2.fastq.gz \
-S ${i%_R1.fastq.gz}.sam
done
