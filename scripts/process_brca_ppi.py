import pandas as pd

brca_ppi = pd.read_csv('breast_cancer_ppi.txt', sep='\t', header=None)
brca_info = pd.read_csv('9606.protein.info.v12.0.txt', sep='\t')

# Drop the third column if it exists
if brca_ppi.shape[1] > 2:
    brca_ppi = brca_ppi.iloc[:, :2]

brca_ppi.columns = ['protein1', 'protein2']

# Filter to protein_id in brca_info only in protein1 and protein2 in brca_ppi
edgelist_proteins = set(brca_ppi['protein1']).union(set(brca_ppi['protein2']))
brca_info_filtered = brca_info[brca_info['protein_id'].isin(edgelist_proteins)]
brca_info_filtered.to_csv('brca_protein_info.csv', index=False)
brca_ppi.to_csv('brca_ppi_edgelist.csv', index=False)
