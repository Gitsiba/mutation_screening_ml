"""Download a specified ENCODE file. Set ENCODE_FILE_ID to a real accession first."""
import os
import re
import urllib.request
from pathlib import Path

file_id = os.environ.get("ENCODE_FILE_ID", "").strip()
if not re.fullmatch(r"ENCFF[0-9A-Z]+", file_id):
    raise SystemExit(
        "Set ENCODE_FILE_ID to a verified ENCODE file accession, e.g. "
        "'ENCFF...' from the ENCODE portal. The repository intentionally "
        "does not ship a placeholder accession."
    )

output_dir = Path(os.environ.get("OUTPUT_DIR", "data_download"))
output_dir.mkdir(parents=True, exist_ok=True)
output_path = output_dir / f"{file_id}.fastq.gz"
url = f"https://www.encodeproject.org/files/{file_id}/@@download/{file_id}.fastq.gz"
request = urllib.request.Request(url, headers={"User-Agent": "mutation-screening-ml/0.1"})
print(f"Downloading {file_id} from ENCODE to {output_path} ...")
try:
    with urllib.request.urlopen(request, timeout=60) as response, output_path.open("wb") as output:
        if response.status != 200:
            raise RuntimeError(f"ENCODE returned HTTP {response.status}")
        while chunk := response.read(1024 * 1024):
            output.write(chunk)
except Exception:
    output_path.unlink(missing_ok=True)
    raise
if output_path.stat().st_size == 0:
    output_path.unlink(missing_ok=True)
    raise SystemExit("Downloaded file is empty; verify the ENCODE accession and file type.")
print(f"Downloaded {output_path} ({output_path.stat().st_size:,} bytes)")
