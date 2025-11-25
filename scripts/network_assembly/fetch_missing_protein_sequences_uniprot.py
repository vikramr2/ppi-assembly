import requests
import csv
import json
import time
from tqdm import tqdm
import pandas as pd
from typing import Optional

def get_protein_sequence(gene_name: str, organism: str = "human") -> Optional[str]:
    """
    Retrieve protein sequence from UniProt using a gene name.
    
    Parameters:
    -----------
    gene_name : str
        The gene name (e.g., 'INPP4A', 'TP53')
    organism : str, optional
        The organism name (default: 'human'). Can use common names like 'human', 'mouse'
        or taxonomy IDs like '9606' for human
    
    Returns:
    --------
    str or None
        The protein sequence in FASTA format, or None if not found
    
    Example:
    --------
    >>> sequence = get_protein_sequence('INPP4A')
    >>> print(sequence[:100])  # Print first 100 characters
    """
    
    # Map common organism names to taxonomy IDs
    organism_map = {
        'human': '9606',
        'mouse': '10090',
        'rat': '10116',
        'yeast': '559292',
        'fly': '7227',
        'worm': '6239'
    }
    
    # Get taxonomy ID
    tax_id = organism_map.get(organism.lower(), organism)
    
    # UniProt REST API endpoint
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    
    # Query parameters - search for gene name and organism
    params = {
        'query': f'(gene:{gene_name}) AND (organism_id:{tax_id})',
        'format': 'fasta',
        'size': 1  # Get only the top result
    }
    
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        
        if response.text.strip():
            return response.text
        else:
            print(f"No protein sequence found for gene '{gene_name}' in {organism}")
            return None
            
    except requests.exceptions.RequestException as e:
        print(f"Error retrieving data from UniProt: {e}")
        return None

def get_protein_sequence_only(gene_name: str, organism: str = "human") -> Optional[str]:
    """
    Retrieve only the amino acid sequence (without FASTA header).
    
    Parameters:
    -----------
    gene_name : str
        The gene name (e.g., 'INPP4A', 'TP53')
    organism : str, optional
        The organism name (default: 'human')
    
    Returns:
    --------
    str or None
        The protein sequence as a single string, or None if not found
    """
    fasta_result = get_protein_sequence(gene_name, organism)
    
    if fasta_result:
        # Split by newlines and skip the header line (starts with '>')
        lines = fasta_result.strip().split('\n')
        sequence = ''.join(line for line in lines if not line.startswith('>'))
        return sequence
    
    return None

current_sequence_mapping = json.load(open('../processed_data/node_to_sequence.json', 'r'))
print(f"Loaded existing mapping with {len(current_sequence_mapping)} entries.")

protein_info_df = pd.read_csv('../processed_data/brca_protein_info.csv')
print(f"Loaded protein info with {len(protein_info_df)} entries.")
print(protein_info_df[['node_id', 'protein_id', 'preferred_name']].head())

# Create a dict mapping node_id to preferred_name (gene name)
protein_info_dict = protein_info_df.set_index('node_id')['preferred_name'].to_dict()
missing_proteins = [node_id for node_id in protein_info_dict if str(node_id) not in current_sequence_mapping]
missing_proteins_dict = {node_id: protein_info_dict[node_id] for node_id in missing_proteins}

# Fetch missing sequences
print(f"Fetching sequences for {len(missing_proteins)} missing proteins using STRING API (gene names)...")
new_mappings = {}
failed_proteins = []

for node_id in tqdm(missing_proteins, desc="Fetching missing sequences"):
    gene_name = missing_proteins_dict[node_id]

    # Fetch sequence
    sequence = get_protein_sequence_only(gene_name)

    if sequence:
        new_mappings[str(node_id)] = sequence
    else:
        failed_proteins.append((node_id, gene_name))

    # Be nice to the API - STRING is generally permissive but let's be respectful
    time.sleep(0.3)

# Update and save mapping to JSON
with open('../processed_data/node_to_sequence_missing_uniprot.json', 'w') as f:
    json.dump(new_mappings, f, indent=2)

print(f"\n=== SUMMARY ===")
print(f"Successfully fetched: {len(new_mappings)}/{len(missing_proteins)}")
print(f"Failed: {len(failed_proteins)}/{len(missing_proteins)}")
print(f"\nNew mapping saved to: ../processed_data/node_to_sequence_missing_uniprot.json")

# Optionally print some failed proteins for debugging
if failed_proteins:
    print(f"\nFirst 10 failed proteins:")
    for node_id, gene_name in failed_proteins[:10]:
        print(f"  Node {node_id}: {gene_name}")
