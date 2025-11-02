# mutation_screening_ml

End-to-end RNA-seq pipeline for cancer mutation screening using ENCODE data, variant calling, and visualization.

## Table of Contents
- [Overview](#overview)
- [Pipeline Steps](#pipeline-steps)
- [Tools Used](#tools-used)
- [How to Run](#how-to-run)
- [Sample Outputs](#sample-outputs)
- [License](#license)

## Overview
This project demonstrates a full workflow for identifying mutations from RNA-seq data. It starts with downloading data from ENCODE, followed by quality control, alignment, variant calling, and visualization of mutation profiles.

## Pipeline Steps
1. Download RNA-seq data from ENCODE
2. Perform quality control and trimming
3. Align reads to reference genome
4. Call variants using bcftools
5. Filter and annotate mutations
6. Generate plots for mutation distribution and impact

## Tools Used
- Python, Bash
- FastQC, fastp
- STAR aligner
- samtools, bcftools
- matplotlib, seaborn

## How to Run

Clone the repo:
```bash
git clone https://github.com/Gitsiba/mutation_screening_ml.git
cd mutation_screening_ml

