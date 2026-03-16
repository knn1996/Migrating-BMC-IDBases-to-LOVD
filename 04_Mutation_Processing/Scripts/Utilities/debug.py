import requests, json

resp = requests.get("https://rest.uniprot.org/uniprotkb/P00813.json", timeout=30)
data = resp.json()

# Print all top-level keys so you can see what's actually there
print("Top-level keys:", list(data.keys()))

# Correct field name is uniProtKBCrossReferences, not dbReferences
embl_refs = [
    ref for ref in data.get("uniProtKBCrossReferences", [])
    if ref.get("database") == "EMBL"
]
print(f"\nTotal EMBL cross-references: {len(embl_refs)}")
for ref in embl_refs[:5]:
    print(json.dumps(ref, indent=2))