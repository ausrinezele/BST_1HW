#!/bin/sh

# Download Genome reference FASTA
wget "https://ftp.ensembl.org/pub/release-67/fasta/mus_musculusls/dna/Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz" 
gzip -d ../../references/Mus_musculus.NCBIM37.67.dna.toplevel.fa.gz

# Download transcriptome reference FASTA.gz
wget "https://ftp.ensembl.org/pub/release-67/fasta/mus_musculus/cdna/Mus_musculus.NCBIM37.67.cdna.all.fa.gz" 
gzip -d ../../references/Mus_musculus.NCBIM37.67.cdna.all.fa.gz

# Download gtf file
wget "https://ftp.ensembl.org/pub/release-67/gtf/mus_musculus/Mus_musculus.NCBIM37.67.gtf.gz"
gzip -d ../../references/Mus_musculus.NCBIM37.67.gtf.gz

# Download raw FASTQ files
prefetch ../../inputs/SRR8985047 ../../inputs/SRR8985048 ../../inputs/SRR8985051 ../../inputs/SRR8985052
fastq-dump --gzip ../../inputs/SRR8985047 ../../inputs/SRR8985048 ../../inputs/SRR8985051 ../../inputs/SRR8985052

 