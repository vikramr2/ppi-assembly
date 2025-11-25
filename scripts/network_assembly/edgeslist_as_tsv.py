import pandas as pd

edgelist_file = '../../processed_data/brca_ppi_edgelist_cleaned.csv'
output_tsv_file = '../../processed_data/brca_ppi_edgelist_cleaned.tsv'

el_df = pd.read_csv(edgelist_file)
el_df.to_csv(output_tsv_file, sep='\t', index=False, header=False)
print(f"Converted {edgelist_file} to TSV format at {output_tsv_file}.")