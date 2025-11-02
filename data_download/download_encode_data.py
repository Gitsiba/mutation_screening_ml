import os
import requests

ENCODE_FASTQ_URL = "https://www.encodeproject.org/files/ENCFF000XXX/@@download/ENCFF000XXX.fastq.gz"
OUTPUT_DIR = "data_download"

os.makedirs(OUTPUT_DIR, exist_ok=True)
output_path = os.path.join(OUTPUT_DIR, "sample.fastq.gz")

print("Downloading ENCODE FASTQ file...")
response = requests.get(ENCODE_FASTQ_URL)
with open(output_path, "wb") as f:
    f.write(response.content)

print(f"Downloaded to {output_path}")

