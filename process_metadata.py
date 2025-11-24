import pandas as pd

enriched_metadata = pd.read_csv('enriched_metadata.csv')
base_metadata = pd.read_csv('metadata.csv')

# Drop gene_type and other_aliases columns from enriched_metadata
enriched_metadata = enriched_metadata.drop(columns=['gene_type', 'other_aliases'])

# Base metadata is sorted in the same order, use its entrez_id and gene_symbol column instead
enriched_metadata['entrez_id'] = base_metadata['entrez_id']
enriched_metadata['gene_symbol'] = base_metadata['gene_name']
enriched_metadata['node_id'] = base_metadata['node_id']


enriched_metadata.to_csv('final_metadata.csv', index=False)
