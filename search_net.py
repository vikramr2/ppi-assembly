import pandas as pd

# Load the extracted file
ppi = pd.read_csv("9606.protein.links.v12.0.txt", sep=' ')

print(f"Total interactions: {len(ppi)}")

# Load protein info
info = pd.read_csv("9606.protein.info.v12.0.txt", sep='\t')

# Check what columns are actually in the info file
print(f"\nInfo file columns: {info.columns.tolist()}")
print(info.head())

# The column is likely 'protein_external_id' or just '#string_protein_id'
# Let's use the correct column name

breast_cancer_genes = [
    # High penetrance
    'BRCA1', 'BRCA2', 'TP53', 'PTEN', 'CDH1', 'STK11',
    # Moderate penetrance
    'ATM', 'CHEK2', 'PALB2', 'BARD1', 'BRIP1', 'RAD51C', 'RAD51D',
    # Hormone receptors
    'ESR1', 'ESR2', 'PGR', 'AR',
    # HER2/EGFR
    'ERBB2', 'EGFR', 'ERBB3',
    # PI3K/AKT
    'PIK3CA', 'PIK3R1', 'AKT1', 'AKT2', 'MTOR',
    # Cell cycle
    'CCND1', 'CDK4', 'CDK6', 'CDKN2A', 'RB1',
    # Other key genes
    'GATA3', 'FOXA1', 'MYC', 'MAP3K1', 'ARID1A'
]

# Map gene names to STRING IDs using the correct column name
bc_proteins = info[info['preferred_name'].isin(breast_cancer_genes)]

# Use the first column (STRING protein ID) - it's likely '#string_protein_id'
protein_id_column = info.columns[0]  # First column is the protein ID
protein_ids = set(bc_proteins[protein_id_column].values)

print(f"\nFound {len(protein_ids)} breast cancer proteins")
print(f"Example IDs: {list(protein_ids)[:3]}")

# Filter for breast cancer interactions (high confidence)
bc_ppi = ppi[
    (ppi['protein1'].isin(protein_ids)) & 
    (ppi['protein2'].isin(protein_ids)) &
    (ppi['combined_score'] > 700)
]

print(f"\nBreast cancer network: {len(bc_ppi)} interactions")

# Save the breast cancer edgelist
bc_ppi.to_csv("breast_cancer_ppi.txt", 
              columns=['protein1', 'protein2', 'combined_score'],
              sep='\t', index=False, header=False)

print("Saved to breast_cancer_ppi.txt")
