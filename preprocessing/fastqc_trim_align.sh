#!/usr/bin/env bash
set -euo pipefail

FASTQ="${FASTQ:-data_download/sample.fastq.gz}"
GENOME_INDEX="${STAR_INDEX:-reference/star_index}"
OUTDIR="${OUTDIR:-preprocessing}"
THREADS="${THREADS:-4}"
TRIMMED="$OUTDIR/trimmed.fastq.gz"
OUT_BAM="$OUTDIR/aligned.bam"

for tool in fastqc fastp STAR samtools; do
  command -v "$tool" >/dev/null 2>&1 || { echo "Required tool not found: $tool" >&2; exit 1; }
done
[[ -s "$FASTQ" ]] || { echo "FASTQ not found or empty: $FASTQ" >&2; exit 1; }
[[ -d "$GENOME_INDEX" ]] || { echo "STAR index directory not found: $GENOME_INDEX" >&2; exit 1; }
mkdir -p "$OUTDIR"

echo "[1/4] FastQC"
fastqc "$FASTQ" -o "$OUTDIR"

echo "[2/4] fastp (single-end input)"
fastp -i "$FASTQ" -o "$TRIMMED" -h "$OUTDIR/fastp.html" -j "$OUTDIR/fastp.json" --thread "$THREADS"

echo "[3/4] STAR alignment"
STAR --runThreadN "$THREADS" --genomeDir "$GENOME_INDEX" \
  --readFilesIn "$TRIMMED" --readFilesCommand zcat \
  --outFileNamePrefix "$OUTDIR/" --outSAMtype BAM SortedByCoordinate

mv "$OUTDIR/Aligned.sortedByCoord.out.bam" "$OUT_BAM"
echo "[4/4] BAM indexing"
samtools index "$OUT_BAM"
echo "Done: $OUT_BAM and $OUT_BAM.bai"
