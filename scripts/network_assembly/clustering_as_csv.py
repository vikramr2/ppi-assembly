import pandas as pd

clustering = pd.read_csv('../../processed_data/brca_ppi_leiden.tsv', sep='\t', header=None, names=['node_id', 'cluster_id'])
clustering.to_csv('../../processed_data/brca_ppi_leiden.csv', index=False)
print(f"Converted '../../processed_data/brca_ppi_leiden.tsv' to CSV format at '../../processed_data/brca_ppi_leiden.csv'.")