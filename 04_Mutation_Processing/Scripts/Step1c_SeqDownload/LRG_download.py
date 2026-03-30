import csv
import os
import requests
from pathlib import Path

_SCRIPT_DIR   = Path(__file__).parent
_MUTATION_DIR = _SCRIPT_DIR.parent.parent

csv_file   = str(_MUTATION_DIR / "Output" / "Step1b_RSG_Mapping" / "LRG.csv")
output_dir = str(_MUTATION_DIR / "DNA sequences" / "LRG FASTA file (Source 3)" / "NG")
os.makedirs(output_dir, exist_ok=True)

base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

with open(csv_file, newline='') as f:
    reader = csv.DictReader(f)

    for row in reader:
        name = row["name"]
        rsg  = row["RSG"]

        fasta_path = os.path.join(output_dir, f"{name}.fasta")
        if os.path.exists(fasta_path):
            print(f"Skipping {name} ({rsg}) — already exists")
            continue

        url = f"{base_url}?db=nuccore&id={rsg}&rettype=fasta&retmode=text"
        print(f"Downloading {name} ({rsg})...")

        response = requests.get(url)

        if response.status_code == 200 and response.text.startswith(">"):
            with open(fasta_path, "w") as out:
                out.write(response.text)
        else:
            print(f"Failed to download {rsg}")
