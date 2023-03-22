#!/bin/bash

threads=6

# Path to reference genome FASTA file
ref_genome="../../references/Mus_musculus.NCBIM37.67.dna.toplevel.fa"
ref_genome_gtf="../../references/Mus_musculus.NCBIM37.67.gtf"
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
for i in ../../outputs/trimmed/*_1_val_1.fq.gz 
do 
    R1=${i}
    R2="../../outputs/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_2_val_2.fq.gz"
    fastqc -t ${threads} ${R1} ${R2} -o ../../outputs/trimmed

done

# Create MultiQC plots for raw and processed data
multiqc ../../outputs -o ../../outputs/multiq_processed

# Mapped each sample to the reference genome. Removed unmapped reads as well as reads that are mapped incorrectly. Mapping results in BAM format
# Deduplicate files and index BAM files.
for i in ../../outputs/trimmed/*_1_val_1.fq.gz
do
    R1=${i}
    R2="../../outputs/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_2_val_2.fq.gz"
    SAM="../../outputs/trimmed/"$(basename ${i} _1_val_1.fq.gz)".sam"
    BAM="../../outputs/trimmed/"$(basename ${i} _1_val_1.fq.gz)".bam"
    BAM_S="../../outputs/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_sorted.bam"
    BAM_Fil="../../outputs/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_filtered.bam"
    BAM_F="../../outputs/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_fixmate.bam"
    BAM_P="../../outputs/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_positionsort.bam"
    BAM_M="../../outputs/trimmed/"$(basename ${i} _1_val_1.fq.gz)"_markdup.bam"
    hisat2 -p ${threads} --dta -x ${ref_genome} -1 ${R1} -2 ${R2} -S ${SAM}
    samtools view -@ 6 -bS ${SAM} -o ${BAM}
    samtools view -@ 6 -F 260 -bS ${BAM} > ${BAM_Fil}
    samtools sort -n -o ${BAM_S} ${BAM_Fil}
    samtools fixmate -m ${BAM_S} ${BAM_F}
    samtools sort -o ${BAM_P} ${BAM_F}
    samtools markdup -r ${BAM_P} ${BAM_M}
    samtools index ${BAM_M}
done

# Stringtie quantification on BAM files.
for i in ../../outputs/trimmed/*_markdup.bam 
do
    GTF="../../outputs/stringtie/"$(basename ${i} _markdup.bam)".gtf"
    BAM_M="../../outputs/trimmed/"$(basename "${i}" .markdup.bam)
    prefix=$(basename ${i} _markdup.bam)
    stringtie -p ${threads} -G ${ref_genome_gtf} -o ${GTF} ${BAM_M}
done

# Create a correlation diagram as well as a PCA plot for data
multiBamSummary bins --outFileName ../../results/mapped.npz --binSize 1000 -p ${threads} --outRawCounts ../../results/raw_counts.tsv -b ../../outputs/trimmed/*_markdup.bam

plotCorrelation -in ../../results/mapped.npz -c pearson -p heatmap -o ../../results/mapped_data_heatmap.pdf
plotPCA -in ../../results/mapped.npz -o ../../results/mapped_data_PCA.pdf



