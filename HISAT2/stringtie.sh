#!/bin/bash
#SBATCH --job-name="trap"
#SBATCH -N 1
#SBATCH --mem=16G
#SBATCH -p medicine_q
#SBATCH --mail-type="ALL"
#SBATCH -t 148:48:15
#SBATCH -J stringtie_ref_guided

for i in *.bam; do

/gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Tools/stringtie-2.1.4.Linux_x86_64/s$
--rf -p 8 \
-G /gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Reference/genome.gtf \
-l ${i%.bam} \
-e \
-B \
-o ${i%.bam}/transcripts.gtf \ 
-A ${i%.bam}/gene_abundances.tsv \
$i
done
