#!/bin/bash
#SBATCH --job-name="trap"
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH -p medicine_q
#SBATCH --mail-type="ALL"
#SBATCH -t 48:48:15
#SBATCH -J MergeBamfiles
java -Xmx2g -jar /gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Tools/picard.jar \
MergeSamFiles OUTPUT=HDF1_UT.bam \
INPUT=HDF1_UT1.bam \
INPUT=HDF1_UT2.bam \
INPUT=HDF1_UT3.bam
