#!/bin/sh

# Calculate (and print to the screen) the number of sequences in your reference genome
grep -c "^>" ../../references/Mus_musculus.NCBIM37.67.cdna.all.fa

# Calculate (and print to the screen) the number of reads in each sample
gzip -cd ../../inputs/SRR8985047.fastq.gz | echo $((`wc -l`/4))
gzip -cd ../../inputs/SRR8985048.fastq.gz | echo $((`wc -l`/4))
gzip -cd ../../inputs/SRR8985051.fastq.gz | echo $((`wc -l`/4))
gzip -cd ../../inputs/SRR8985052.fastq.gz | echo $((`wc -l`/4))

# Calculate the number of protein-coding genes in your genome
awk '$2 == "protein_coding"' ../../references/Mus_musculus.NCBIM37.67.gtf.gz | cut -f 9 | cut -f 2 -d ";" | sort -u | wc -l
