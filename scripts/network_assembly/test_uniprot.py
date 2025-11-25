import requests
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


# Example usage
if __name__ == "__main__":
    # Test with INPP4A
    print("Fetching INPP4A protein sequence...")
    sequence = get_protein_sequence('INPP4A')
    
    if sequence:
        print("\nFull FASTA format:")
        print(sequence[:200] + "..." if len(sequence) > 200 else sequence)
        
        # Get just the sequence
        seq_only = get_protein_sequence_only('INPP4A')
        print(f"\nSequence length: {len(seq_only)} amino acids")
        print(f"First 50 amino acids: {seq_only[:50]}")
    
    # Try another gene
    print("\n" + "="*50)
    print("\nFetching TP53 protein sequence...")
    tp53_seq = get_protein_sequence_only('TP53')
    if tp53_seq:
        print(f"TP53 sequence length: {len(tp53_seq)} amino acids")
        print(f"First 50 amino acids: {tp53_seq[:50]}")
        