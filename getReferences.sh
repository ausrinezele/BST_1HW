#!/bin/sh

efetch -db SRA -id SRR8985047 -format fasta > SRR8985047.fasta
efetch -db SRA -id SRR8985048 -format fasta > SRR8985048.fasta
efetch -db SRA -id SRR8985051 -format fasta > SRR8985051.fasta
efetch -db SRA -id SRR8985052 -format fasta > SRR8985052.fasta

wget -O SRR8985047.gtf "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=SRR8985047"
wget -O SRR8985048.gtf "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=SRR8985048"
wget -O SRR8985051.gtf "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=SRR8985052"
wget -O SRR8985052.gtf "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=SRR8985052"

fastqc inputs/SRR8985047.fastq.gz
