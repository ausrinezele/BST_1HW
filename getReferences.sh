#!/bin/sh

# Download Genome reference FASTA
wget "https://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz" 
gzip -d Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz

# Download transcriptome reference FASTA.gz
wget "https://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.67.cdna.all.fa.gz" 
gzip -d Mus_musculus.NCBIM37.67.cdna.all.fa.gz

# Download gtf file
wget "https://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/dna/Mus_musculus.NCBIM37.67.gtf.gz"
gzip -d Mus_musculus.NCBIM37.67.gtf.gz

fastq-dump --outdir /inputs --gzip --skip-technical --readids --dumpbase --split-files SRR8985047 SRR8985048 SRR8985051 SRR8985052

 