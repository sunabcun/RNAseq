#!/bin/bash
#SBATCH --job-name="trap"
#SBATCH -n 5
#SBATCH -p medicine_q
#SBATCH --mail-type="ALL"
#SBATCH -t 48:48:15
#SBATCH -J FASTQC
module load fastqc
fastqc -t 6 *.fastq.gz

mv *.html ~/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples2020_Outputs/QC/
mv *.zip
~/2020_Outputs/Akash_Gunjan_05-19-2020_Yuna-samples2020_Outputs/QC/
