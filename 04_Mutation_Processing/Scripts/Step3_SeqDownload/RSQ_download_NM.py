import csv
import os
import requests
from pathlib import Path

_SCRIPT_DIR   = Path(__file__).parent
_MUTATION_DIR = _SCRIPT_DIR.parent.parent

csv_file   = str(_MUTATION_DIR / "Output" / "LRG_with_NM.csv")
output_dir = str(_MUTATION_DIR / "DNA sequences" / "LRG FASTA file (Source 3)" / "NM")
os.makedirs(output_dir, exist_ok=True)

base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

with open(csv_file, newline='') as f:
    reader = csv.DictReader(f)

    for row in reader:
        name = row["name"]
        nm   = row["NM"]

        url = f"{base_url}?db=nuccore&id={nm}&rettype=fasta&retmode=text"

        print(f"Downloading {name} ({nm})...")

        response = requests.get(url)

        if response.status_code == 200 and response.text.startswith(">"):
            fasta_path = os.path.join(output_dir, f"{name}_{nm}.fasta")
            with open(fasta_path, "w") as out:
                out.write(response.text)
            print(f"  Saved {name}_{nm}.fasta")
        else:
            print(f"  Failed to download {nm}")
