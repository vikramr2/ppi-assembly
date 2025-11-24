import pandas as pd

brca_info = pd.read_csv('brca_protein_info.csv')
brca_edges = pd.read_csv('brca_ppi_edgelist.csv')

protein_list = dict(enumerate(brca_info['protein_id'].tolist()))
brca_info['node_id'] = brca_info['protein_id'].map({v: k for k, v in protein_list.items()})

brca_edges['source'] = brca_edges['protein1'].map({v: k for k, v in protein_list.items()})
brca_edges['target'] = brca_edges['protein2'].map({v: k for k, v in protein_list.items()})

brca_info.to_csv('brca_protein_info.csv', index=False)
brca_edges[['source', 'target']].to_csv('brca_ppi_edgelist.csv', index=False)