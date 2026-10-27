"""Call candidate variants from one aligned RNA-seq BAM using bcftools.

This is a research prototype. RNA-seq calls require careful artifact review and
are not equivalent to validated somatic mutation calls.
"""
import os
import shutil
import subprocess
from pathlib import Path

BAM = Path(os.environ.get("BAM", "preprocessing/aligned.bam"))
REF = Path(os.environ.get("REFERENCE_FASTA", "reference/genome.fa"))
VCF = Path(os.environ.get("VCF_OUTPUT", "mutation_screening/variants.vcf"))
THREADS = os.environ.get("THREADS", "4")

for tool in ("samtools", "bcftools"):
    if shutil.which(tool) is None:
        raise SystemExit(f"Required executable not found on PATH: {tool}")
for path in (BAM, REF):
    if not path.is_file():
        raise SystemExit(f"Required input file not found: {path}")

VCF.parent.mkdir(parents=True, exist_ok=True)
print(f"Indexing BAM: {BAM}")
subprocess.run(["samtools", "index", str(BAM)], check=True)
print("Running bcftools mpileup and variant calling...")
with subprocess.Popen(
    ["bcftools", "mpileup", "--threads", THREADS, "-Ou", "-f", str(REF), str(BAM)],
    stdout=subprocess.PIPE,
) as mpileup:
    try:
        with VCF.open("w") as output:
            subprocess.run(
                ["bcftools", "call", "--threads", THREADS, "-mv", "-Ov"],
                stdin=mpileup.stdout,
                stdout=output,
                check=True,
            )
    finally:
        if mpileup.stdout:
            mpileup.stdout.close()
    return_code = mpileup.wait()
    if return_code:
        VCF.unlink(missing_ok=True)
        raise subprocess.CalledProcessError(return_code, mpileup.args)
print(f"Wrote candidate VCF: {VCF}")
print("Reminder: review RNA-seq-specific artifacts; this output is not clinically validated.")
