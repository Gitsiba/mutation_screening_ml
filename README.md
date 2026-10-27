# Mutation Screening Workbench

Research-oriented RNA-seq variant-calling prototype with a browser-based VCF explorer. Built with Python, Bash, FastQC, fastp, STAR, samtools, bcftools, and a lightweight client-side dashboard.

![Mutation ML dashboard preview](docs/mutation_ml_preview.svg)

> **Status:** Research prototype. This repository does **not** currently contain a trained machine-learning model, validated somatic mutation classifier, or clinical interpretation system. The preview image is conceptual; dashboard metrics and charts are calculated only when a user loads a real VCF into the browser app.

[**Open the browser VCF workbench**](https://gitsiba.github.io/mutation_screening_ml/) · [Source code](https://github.com/Gitsiba/mutation_screening_ml)

## Why this project

The project demonstrates an auditable path from sequencing reads to candidate variant records and exploratory summaries. The browser app parses a selected plain-text VCF locally, calculates row counts and chromosome distributions from the loaded file, supports search and FILTER=PASS filtering, and exports the filtered rows to CSV. The selected file is not uploaded to a server.

## Current capabilities

- Configurable ENCODE FASTQ download using a real file accession supplied by the user.
- Single-end read QC, trimming, STAR alignment, and BAM indexing script.
- bcftools mpileup/call wrapper that uses a real subprocess pipe and checks failures.
- Browser-based VCF explorer with local parsing, searchable records, FILTER filtering, chromosome counts, and CSV export.
- Research-use warnings and documented input/reference assumptions.

## Workflow overview

```text
ENCODE FASTQ accession
        |
        v
FastQC -> fastp -> STAR -> sorted BAM -> samtools index
                                          |
Reference FASTA + index ------------------+
                                          v
                               bcftools mpileup -> bcftools call
                                          |
                                          v
                                      variants.vcf
                                          |
                                          v
                              Browser VCF exploration
```

## Quick start

### 1. Download a real ENCODE file

Find a suitable FASTQ file in the [ENCODE portal](https://www.encodeproject.org/). Confirm organism, assay, biosample, read layout, and accession before downloading. The example intentionally uses no fabricated accession.

```bash
python data_download/download_encode_data.py
```

Set the accession first (replace the example text with a real file ID):

```bash
export ENCODE_FILE_ID=ENCFF_REPLACE_WITH_REAL_ACCESSION
python data_download/download_encode_data.py
```

The downloader accepts an ENCODE file accession matching `ENCFF...`; it does not establish that the selected file is biologically appropriate for mutation calling.

### 2. Prepare references

Create a STAR genome index matching the selected sample and reference build, and place it at `reference/star_index/` (or set `STAR_INDEX`). Place the matching reference FASTA at `reference/genome.fa` and index it with `samtools faidx reference/genome.fa`. Reference genome, annotation, and build must match the sample and analysis plan.

### 3. Run single-end QC and alignment

The included shell script is configured for **single-end** FASTQ input. It is not a paired-end workflow.

```bash
FASTQ=data_download/REAL_FILE_ID.fastq.gz STAR_INDEX=reference/star_index THREADS=4 bash preprocessing/fastqc_trim_align.sh
```

Required tools: FastQC, fastp, STAR, and samtools. The STAR index and reference must already exist.

### 4. Call candidate variants

```bash
REFERENCE_FASTA=reference/genome.fa BAM=preprocessing/aligned.bam \
  VCF_OUTPUT=mutation_screening/variants.vcf python mutation_screening/variant_calling_pipeline.py
```

Required tools: bcftools and samtools. Review the exact caller parameters and output for the research question before drawing conclusions.

### 5. Explore the VCF in a browser

Open the [Mutation Screening Workbench](https://gitsiba.github.io/mutation_screening_ml/), choose a plain-text `.vcf` file, search or filter records, inspect the chromosome count chart, and export the visible filtered data. Compressed `.vcf.gz` is not supported by the browser demo; decompress it locally first. Large files can use substantial browser memory.

## Repository layout

```text
data_download/download_encode_data.py        Configurable ENCODE file downloader
preprocessing/fastqc_trim_align.sh            Single-end QC, trimming, alignment, BAM indexing
mutation_screening/variant_calling_pipeline.py bcftools candidate variant calling
visualization/mutation_plots.ipynb            Notebook for exploratory variant counts
docs/index.html                               Client-side VCF workbench
docs/mutation_ml_preview.svg                  Conceptual interface preview
```

## Interpretation and limitations

- **RNA-seq variant calls are candidate calls, not confirmed mutations.** Mapping artifacts, RNA editing, splice junctions, allele-specific expression, coverage, and transcript expression can affect results.
- A single RNA-seq sample and a VCF alone do not establish somatic status; matched normal DNA, suitable controls, validated filters, and orthogonal evidence may be required for the intended study.
- The pipeline is not a complete production workflow: no pinned environment/container, automated sample metadata validation, paired-end orchestration, comprehensive annotation, benchmarking, or end-to-end test dataset is included yet.
- No trained ML model, model weights, training labels, evaluation metrics, or external validation are included. The dashboard does not make predictions.
- This is not a clinical device or diagnostic tool. Do not use it for patient-care or treatment decisions.

## Reproducibility and roadmap

Next engineering steps: add a pinned environment, test fixtures and automated tests; configure sample metadata and paired-end support if needed; document reference-build provenance; add validated annotation; and only then evaluate a clearly defined ML task using curated labels, train/test separation, class-imbalance handling, calibration, and independent validation.

## License

No license is currently specified. Add an explicit license before inviting reuse or external contributions.
