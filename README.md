# RNAseq

# 1. Trimming the fastq files and get QC data.
  - TrimGalore: https://github.com/FelixKrueger/TrimGalore
  - To run this, need to have cutadapt and fastq
  - cutadapt: https://cutadapt.readthedocs.io/en/stable/installation.html
          python3 -m pip install --user --upgrade cutadapt
  - fastqc
    - wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
    - unzip fastqc_v0.11.9.zip
    - cd FastQC/
    - chmod 755 fastqc
    - ./fastqc --help
  - multiqc: A tool for assembling QC reports is a python package
    - pip3 install --user multiqc
    - export PATH=/home/ubuntu/.local/bin:$PATH
    - python3 -m multiqc --help

# 2-1. Align (STAR)
  - Index genome for STAR
    - wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
    - gunzip GRCh38.primary_assembly.genome.fa.gz
    - FASTA="../GRCh38.primary_assembly.genome.fa"
    - wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/gencode.v34.primary_assembly.annotation.gtf.gz
    - gunzip gencode.v34.primary_assembly.annotation.gtf.gz
    - GTF="../gencode.v34.primary_assembly.annotation.gtf"
  - Generate genome indices:index.sh
  - Mapping to the reference genome: map.sh
 
 # 2-2. Align (SALMON)
 - Quantify abundance and effective transcript lengths
 - Download transcripts info: gencode.v37.transcripts.fa.gz
 - Build index: index.sh
 - Mapping 
# 3. Post alignment QC
  - flagstat: generate alignment metrics
    - flagstat.sh
  - fastqc + multiqc: multiqc.sh
