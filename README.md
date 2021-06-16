# RNAseq

# 1. Trimming the fastq files and get QC data.
  - TrimGalore: https://github.com/FelixKrueger/TrimGalore
  - To run this, need to have cutadapt and fastq
  - cutadapt: https://cutadapt.readthedocs.io/en/stable/installation.html
          python3 -m pip install --user --upgrade cutadapt

# 2. Align (STAR)
  - Index genome for STAR
    - wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
    - gunzip GRCh38.primary_assembly.genome.fa.gz
    - FASTA="../GRCh38.primary_assembly.genome.fa"
    - wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz
    - gunzip gencode.v34.primary_assembly.annotation.gtf.gz
    - GTF="../gencode.v34.primary_assembly.annotation.gtf"
  - Generate genome indices:index.sh
  - Mapping to the reference genome: map.sh
