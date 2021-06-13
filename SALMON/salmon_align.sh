#!/bin/bash
#SBATCH --job-name="trap"
#SBATCH -N 1
#SBATCH --mem=24G
#SBATCH -p medicine_q
#SBATCH --mail-type="ALL"
#SBATCH -t 148:48:15
#SBATCH -J salmon
for i in *_R1.fq.gz; do
/gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Tools/salmon/bin/salmon quant \
--threads 8 \
--index /gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Reference/Salmon/GRCh38_salmon_index \
--libType A \
--validateMappings \
--geneMap /gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Reference/Salmon/gencode.v37.annotation.gtf \
--output /gpfs/research/medicine/sequencer/NovaSeq/Outputs_fastq/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples/HISAT2/Salmon/${i%_R1.fq.gz} \
-1 $i \
-2 ${i%_R1.fq.gz}_R2.fq.gz
done
