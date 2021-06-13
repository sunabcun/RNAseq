#!/bin/bash
#SBATCH --job-name="trap"
#SBATCH -N 6
#SBATCH --mem=16G
#SBATCH -p medicine_q
#SBATCH --mail-type="ALL"
#SBATCH -t 148:48:15
#SBATCH -J SAMtoBAM
for i in *.sam; do
/gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Tools/samtools-1.11/samtools sort -@ 8 \
-o ${i%.sam}.bam \
$i
done
