import requests

def get_uniprot_sequence(ensp_id):
    """Fetch sequence from UniProt using Ensembl ID"""
    # Remove "9606." prefix if present
    ensp_id = ensp_id.split('.')[-1]
    
    # Query UniProt for this Ensembl ID
    url = f"https://rest.uniprot.org/uniprotkb/search"
    params = {
        'query': f'xref:ensembl-{ensp_id}',
        'format': 'fasta'
    }
    
    response = requests.get(url, params=params)
    
    if response.ok:
        # Parse FASTA format
        lines = response.text.strip().split('\n')
        sequence = ''.join(lines[1:])  # Skip header line
        return sequence
    return None

# Usage
seq = get_uniprot_sequence("9606.ENSP00000000233")
print(seq)
