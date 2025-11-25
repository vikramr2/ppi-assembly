import requests
import csv
import json
import time
from tqdm import tqdm
import pandas as pd

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
    
current_sequence_mapping = json.load(open('../processed_data/node_to_sequence.json', 'r'))
print(f"Loaded existing mapping with {len(current_sequence_mapping)} entries.")

protein_info_df = pd.read_csv('../processed_data/brca_protein_info.csv')
print(f"Loaded protein info with {len(protein_info_df)} entries.")
print(protein_info_df[['node_id', 'protein_id']].head())
protein_info_dict = protein_info_df.set_index('node_id')['protein_id'].to_dict()
missing_proteins = [node_id for node_id in protein_info_dict if str(node_id) not in current_sequence_mapping]
missing_proteins_dict = {node_id: protein_info_dict[node_id] for node_id in missing_proteins}

# Fetch missing sequences
print(f"Fetching sequences for {len(missing_proteins)} missing proteins using Ensembl API...")
new_mappings = {}
failed_proteins = []

for node_id in tqdm(missing_proteins, desc="Fetching missing sequences"):
    protein_id = missing_proteins_dict[node_id]

    # Fetch sequence
    sequence = get_ensembl_sequence(protein_id)

    if sequence:
        new_mappings[str(node_id)] = sequence
    else:
        failed_proteins.append((node_id, protein_id))

    # Be nice to the API - Ensembl allows up to 15 requests per second
    time.sleep(0.5)

# Update and save mapping to JSON
with open('../processed_data/node_to_sequence_missing.json', 'w') as f:
    json.dump(new_mappings, f, indent=2)
    
print(f"\n=== SUMMARY ===")
print(f"Successfully fetched: {len(new_mappings)}/{len(missing_proteins)}")
print(f"Failed: {len(failed_proteins)}/{len(missing_proteins)}")
print(f"\nNew mapping saved to: ../processed_data/node_to_sequence_missing.json")
