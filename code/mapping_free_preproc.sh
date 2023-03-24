#!/bin/bash

threads=6
# Reference transcriptome preparation

transcriptome_dir="../references/Mus_musculus.NCBIM37.67.cdna.all.fa"
transcriptome_index_check="../references/mm9_cDNA/ref_indexing.log"
transcriptome_index_dir="../references/mm9_cDNA"
if [ -f "$transcriptome_index_check" ]
then
    echo "Reference transcriptome index already exists"
else 
    salmon index -t ${transcriptome_dir} -i ${transcriptome_index_dir}
    echo "Reference transcriptome index created"
fi

# Run FASTQC analysis on each of your FASTQC files
for i in ../inputs/*_1.fastq.gz; 
do 
    R1=${i}
    R2="../inputs/"$(basename ${i} _1.fastq.gz)"_2.fastq.gz"
    fastqc -t ${threads} ${R1} ${R2} -o ../outputs/without_mapping

done

# Generate MULTIQC report for results
multiqc ../outputs/without_mapping -o ../outputs/without_mapping/multiq

# Run standard FASTQ trimming
trim_galore -j ${threads} --length 20 -o ../outputs/without_mapping/trimmed --paired  ../inputs/SRR8985047_1.fastq.gz  ../inputs/SRR8985047_2.fastq.gz 
trim_galore -j ${threads} --length 20 -o ../outputs/without_mapping/trimmed --paired  ../inputs/SRR8985048_1.fastq.gz  ../inputs/SRR8985048_2.fastq.gz
trim_galore -j ${threads} --length 20 -o ../outputs/without_mapping/trimmed --paired  ../inputs/SRR8985051_1.fastq.gz  ../inputs/SRR8985051_2.fastq.gz
trim_galore -j ${threads} --length 20 -o ../outputs/without_mapping/trimmed --paired  ../inputs/SRR8985052_1.fastq.gz  ../inputs/SRR8985052_2.fastq.gz

# Rerun FASTQC on newly created/cleaned FASTQ files
for i in ../outputs/without_mapping/trimmed/*_1_val_1.fq.gz 
do 
    R1=${i}
    R2="../outputs/without_mapping/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_2_val_2.fq.gz"
    fastqc -t ${threads} ${R1} ${R2} -o ../outputs/without_mapping/trimmed

done

# Create MultiQC plots for raw and processed data
multiqc ../outputs/without_mapping -o ../outputs/without_mapping/multiq_processed

# Quantification
for i in ../outputs/without_mapping/trimmed/*_1_val_1.fq.gz 
do
    R1=${i}
    R2="../outputs/without_mapping/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_2_val_2.fq.gz"
    QUANT="../outputs/without_mapping/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_transcripts_quant"
    salmon quant -i ${transcriptome_index_dir} -l IU -1 ${R1} -2 ${R2} --validateMappings -o ${QUANT}
done

