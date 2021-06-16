#!/bin/bash
#SBATCH --job-name="trap"
#SBATCH -N 3
#SBATCH --mem=48G
#SBATCH -p medicine_q
#SBATCH --mail-type="ALL"
#SBATCH -t 48:48:15
#SBATCH -J Indexingbam
for i in *.bam; do
/gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Tools/samtools-1.11/samtools index $i
done
