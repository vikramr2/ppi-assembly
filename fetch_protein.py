import requests
import csv
import json
import time
from tqdm import tqdm

def get_ensembl_sequence(ensp_id):
    """Fetch sequence from Ensembl using Ensembl protein ID"""
    # Remove "9606." prefix if present
    ensp_id = ensp_id.split('.')[-1]

    # Use Ensembl REST API
    url = f"https://rest.ensembl.org/sequence/id/{ensp_id}"
    headers = {"Content-Type": "application/json"}

    try:
        response = requests.get(url, headers=headers, timeout=10)

        if response.ok:
            data = response.json()
            return data.get('seq', None)
        return None
    except Exception as e:
        return None

# Read CSV and create mapping
node_to_sequence = {}
failed_proteins = []

with open('brca_protein_info.csv', 'r') as f:
    reader = csv.DictReader(f)
    proteins = list(reader)
    total = len(proteins)

print(f"Processing {total} proteins using Ensembl API...")

for row in tqdm(proteins, desc="Fetching sequences"):
    node_id = int(row['node_id'])
    protein_id = row['protein_id']

    # Fetch sequence
    sequence = get_ensembl_sequence(protein_id)

    if sequence:
        node_to_sequence[node_id] = sequence
    else:
        failed_proteins.append((node_id, protein_id))

    # Be nice to the API - Ensembl allows up to 15 requests per second
    time.sleep(0.07)

# Save mapping to JSON
with open('node_to_sequence.json', 'w') as f:
    json.dump(node_to_sequence, f, indent=2)

print(f"\n=== SUMMARY ===")
print(f"Successfully fetched: {len(node_to_sequence)}/{total}")
print(f"Failed: {len(failed_proteins)}/{total}")
print(f"\nMapping saved to: node_to_sequence.json")

if failed_proteins:
    print(f"\nFailed proteins (first 10):")
    for node_id, protein_id in failed_proteins[:10]:
        print(f"  Node {node_id}: {protein_id}")
    if len(failed_proteins) > 10:
        print(f"  ... and {len(failed_proteins) - 10} more")
