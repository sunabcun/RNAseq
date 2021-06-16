#!/bin/bash
#SBATCH --job-name="trap"
#SBATCH -N 1
#SBATCH --mem=64G
#SBATCH -p medicine_q
#SBATCH --mail-type="ALL"
#SBATCH -t 48:48:15
#SBATCH -J STARalign

for i in *_R1.fastq.gz; do
STAR --runThreadN 8 \
--genomeDir /gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/results/STAR/References \
--outSAMtype BAM SortedByCoordinate \
--quantMode GeneCounts \
--readFilesCommand zcat \
--readFilesIn $i ${i%_R1.fastq.gz}_R2.fastq.gz \
--outFileNamePrefix ${i%_R1.fastq.gz}
done
