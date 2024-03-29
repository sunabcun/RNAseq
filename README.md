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
 - Build index: build_index.sh
 - Align: salmon_align.sh

# 2-3. Align (HISAT2)
- Download pre-built index
  - https://daehwankimlab.github.io/hisat2/download/
- Align: HISAT2align.sh
- SAMtoBAM conversion: sam_bam_conversion.sh
- merge bam files for IGV: merge.sh
- post QC
- stringtie(ref_guided): stringtie.sh
- Ballgown for DE: ballgown.sh
- index bamfiles for HTseq: indexbam.sh
- HTseq count: Run htseq-count on alignments instead to produce raw counts instead of FPKM/TPM values for differential expression analysis
  - HT_seq_count.sh

# 3. Post alignment QC
  - flagstat: generate alignment metrics
    - flagstat.sh
  - fastqc + multiqc: multiqc.sh

# 4. DE analysis/isoform analysis/Gene expression plotting
- Easy to run using Salmon files
- general_DE_gene_exp_plotting.Rmd
- R codes were modified from 3DRNAseq and many previous works
- 3DRNAseq codes were good to modify for isoform analysis and DTU etc.
  - https://www.biorxiv.org/content/10.1101/656686v1

# miRNAseq
- DESEQ2 and plotting: miRNAseq.Rmd

# miRNA + WGBS + RNAseq
- WGBS + RNAseq: 
  - Combine the data with matched gene symbol
  - Remove the rows if they are not matched with higher methylation + low exp or low methylation + high exp (logfc +/- 1)
  - Filter the data by FDR < 0.05 from RNAseq
  - "Low_expressed_genes_in_Keloids.csv""High_expressed_genes_in_Keloids.csv""Low_expressed_genes_name_in_Keloids.csv""High_expressed_genes_name_in_Keloids.csv"
  - Run Cluego to check the overlapping pathways: High_exp_keloids.cluego.cys, Low_exp_keloids.cluego.cys

- WGBS + miRNAseq
  - miRNA has different annotations so had to change the symbol to official names using biomart: miR_name_biomart.txt
  - Modify name to recognize each other between WGBS and miRNAseq data: comparing_difference_wgbs_miRNA_RNAseq.Rmd
  - Remove the rows if they are not matched with higher methylation(low miRNA exp) or low methylation(high miRNA exp). 
  - "Low_expressed_miR_in_Keloids.csv", "High_expressed_miR_in_Keloids.csv", "Low_expressed_miRname_in_Keloids.csv", "High_expressed_miRname_in_Keloids.csv"

- Using High/Low expressed miRNA name in Keloids, search the genes related to the miRNA.
  - The genes are retrieved from miRNA gene data base (miRTarBase release 8.0) using cytargetlinker.
  - https://cytargetlinker.github.io/pages/linksets/mirtarbase
  - High_exp_miRNAgenes.csv, Low_exp_miRNAgenes.csv
  - Merge the genes if they are matched with RNAseq data.
  - Filter the data if highly expressed miRNA related RNAseq genes has lower expression in keloids. (High miRNA -> regulate gene expression(low exp)).
    - Filter HighMir with RNAseq logFC >= 1 (lower expressed genes in keloids (2 times lower in keloids == logFC >= 1))
  - "HighmiRNA_lowgenes_inkeloids.csv", "LowmiRNA_lowgenes_inkeloids.csv"
  - Use Cytoscape to make a figure for miRNA and related gene connections: https://cytargetlinker.github.io/cytoscape-tutorials/protocols/modules/cytargetlinker-tutorial1/#/0/11
  - Use the gene name for checking overlapped GO functions (ClueGO).
 
