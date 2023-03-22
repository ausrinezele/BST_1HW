#!/bin/bash

threads=6
# Reference transcriptome preparation

transcriptome_dir="../../references/Mus_musculus.NCBIM37.67.cdna.all.fa"
transcriptome_dir_index="../../references/Mus_musculus.NCBIM37.67.cdna.all_index"

if [ ! -f "${transcriptome_dir_index}.1.ht2" ] ; then
    hisat2-build -p ${threads} ${transcriptome_dir}  ${transcriptome_dir_index}
    echo "Reference transcriptome index created."
else
    echo "Reference transcriptome index files already exist."
fi

# Run FASTQC analysis on each of your FASTQC files
for i in ../../inputs/*_1.fastq.gz; 
do 
    R1=${i}
    R2="../../inputs/"$(basename ${i} _1.fastq.gz)"_2.fastq.gz"
    fastqc -t ${threads} ${R1} ${R2} -o ../../outputs/without_mapping

done

# Generate MULTIQC report for results
multiqc ../../outputs/without_mapping -o ../../outputs/without_mapping/multiq

# Run standard FASTQ trimming
trim_galore -j ${threads} --length 20 -o ../../outputs/without_mapping/trimmed --paired  ../../inputs/SRR8985047_1.fastq.gz  ../../inputs/SRR8985047_2.fastq.gz
trim_galore -j ${threads} --length 20 -o ../../outputs/without_mapping/trimmed --paired  ../../inputs/SRR8985048_1.fastq.gz  ../../inputs/SRR8985048_2.fastq.gz
trim_galore -j ${threads} --length 20 -o ../../outputs/without_mapping/trimmed --paired  ../../inputs/SRR8985051_1.fastq.gz  ../../inputs/SRR8985051_2.fastq.gz
trim_galore -j ${threads} --length 20 -o ../../outputs/without_mapping/trimmed --paired  ../../inputs/SRR8985052_1.fastq.gz  ../../inputs/SRR8985052_2.fastq.gz

# Rerun FASTQC on newly created/cleaned FASTQ files
for i in ../../outputs/without_mapping/trimmed/*_1_val_1.fq.gz 
do 
    R1=${i}
    R2="../../outputs/without_mapping/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_2_val_2.fq.gz"
    fastqc -t ${threads} ${R1} ${R2} -o ../../outputs/without_mapping/trimmed

done

# Create MultiQC plots for raw and processed data
multiqc ../../outputs/without_mapping -o ../../outputs/without_mapping/multiq_processed