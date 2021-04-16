# RNAseq

1. Trimming the fastq files and get QC data.
  - TrimGalore: https://github.com/FelixKrueger/TrimGalore
  - To run this, need to have cutadapt and fastq
  - cutadapt: https://cutadapt.readthedocs.io/en/stable/installation.html
          python3 -m pip install --user --upgrade cutadapt

2. Align
  - Download prebuilt index: https://daehwankimlab.github.io/hisat2/download/
  - Align with HISAT2
