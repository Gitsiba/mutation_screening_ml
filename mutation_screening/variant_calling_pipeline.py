import subprocess

BAM = "preprocessing/aligned.bam"
VCF = "mutation_screening/variants.vcf"
REF = "reference/genome.fa"

subprocess.run(["mkdir", "-p", "mutation_screening"])

print("Indexing BAM...")
subprocess.run(["samtools", "index", BAM])

print("Calling variants with bcftools...")
subprocess.run([
    "bcftools", "mpileup", "-Ou", "-f", REF, BAM,
    "|", "bcftools", "call", "-mv", "-Ov", "-o", VCF
])
