#!/bin/bash

threads=6

# Path to reference genome FASTA file
ref_genome="../../references/Mus_musculus.NCBIM37.67.dna.toplevel.fa"

# check if HISAT2 index files exist for reference genome
if [ ! -f ${ref_genome}.1.ht2 ]; then
  echo "HISAT2 index files not found. Generating index..."
  hisat2-build -p ${threads} ${ref_genome} ${ref_genome}
  echo "HISAT2 index files generated."
else
  echo "HISAT2 index files already exist."
fi

# Run FASTQC analysis on each of your FASTQC files
for i in ../../inputs/*_1.fastq.gz; 
do 
    R1=${i}
    R2="../../inputs/"$(basename ${i} _1.fastq.gz)"_2.fastq.gz"
    fastqc -t ${threads} ${R1} ${R2} -o ../../outputs

done

# Generate MULTIQC report for results
multiqc ../../outputs -o ../../outputs/multiq

# Run standard FASTQ trimming
trim_galore -j ${threads} --length 20 -o ../../outputs/trimmed --paired  ../../inputs/SRR8985047_1.fastq.gz  ../../inputs/SRR8985047_2.fastq.gz
trim_galore -j ${threads} --length 20 -o ../../outputs/trimmed --paired  ../../inputs/SRR8985048_1.fastq.gz  ../../inputs/SRR8985048_2.fastq.gz
trim_galore -j ${threads} --length 20 -o ../../outputs/trimmed --paired  ../../inputs/SRR8985051_1.fastq.gz  ../../inputs/SRR8985051_2.fastq.gz
trim_galore -j ${threads} --length 20 -o ../../outputs/trimmed --paired  ../../inputs/SRR8985052_1.fastq.gz  ../../inputs/SRR8985052_2.fastq.gz

# Rerun FASTQC on newly created/cleaned FASTQ files
for i in ../../outputs/trimmed/*_1_val_1.fq.gz; 
do 
    R1=${i}
    R2="../../outputs/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_2_val_2.fq.gz"
    fastqc -t ${threads} ${R1} ${R2} -o ../../outputs/trimmed

done

# Create MultiQC plots for raw and processed data
multiqc ../../outputs -o ../../outputs/multiq_processed