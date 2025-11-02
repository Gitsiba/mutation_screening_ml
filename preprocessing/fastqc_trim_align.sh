#!/bin/bash

FASTQ="data_download/sample.fastq.gz"
TRIMMED="preprocessing/trimmed.fastq.gz"
GENOME_INDEX="reference/star_index"
OUT_BAM="preprocessing/aligned.bam"

mkdir -p preprocessing

echo "Running FastQC..."
fastqc $FASTQ -o preprocessing/

echo "Trimming with fastp..."
fastp -i $FASTQ -o $TRIMMED --detect_adapter_for_pe -h preprocessing/fastp.html

echo "Aligning with STAR..."
STAR --runThreadN 4 \
     --genomeDir $GENOME_INDEX \
     --readFilesIn $TRIMMED \
     --readFilesCommand zcat \
     --outFileNamePrefix preprocessing/ \
     --outSAMtype BAM SortedByCoordinate

mv preprocessing/Aligned.sortedByCoord.out.bam $OUT_BAM
