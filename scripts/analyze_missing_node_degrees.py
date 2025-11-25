import json
import pandas as pd
from collections import Counter

# Load the current sequence mapping
current_sequence_mapping = json.load(open('../processed_data/node_to_sequence.json', 'r'))
print(f"Loaded existing mapping with {len(current_sequence_mapping)} entries.")

# Load protein info to get all nodes
protein_info_df = pd.read_csv('../processed_data/brca_protein_info.csv')
print(f"Loaded protein info with {len(protein_info_df)} entries.")

# Find missing nodes
protein_info_dict = protein_info_df.set_index('node_id')['protein_id'].to_dict()
missing_nodes = [node_id for node_id in protein_info_dict if str(node_id) not in current_sequence_mapping]
print(f"Found {len(missing_nodes)} nodes without sequences.")

# Load the edge list
edgelist_df = pd.read_csv('../processed_data/brca_ppi_edgelist_cleaned.csv')
print(f"Loaded edge list with {len(edgelist_df)} edges.")

# Calculate degree for all nodes
degree_counter = Counter()

# Count edges for each node (both as source and target)
for _, row in edgelist_df.iterrows():
    degree_counter[row['source']] += 1
    degree_counter[row['target']] += 1

# Get degrees for missing nodes
missing_node_degrees = []
for node_id in missing_nodes:
    degree = degree_counter.get(node_id, 0)
    protein_id = protein_info_dict[node_id]
    missing_node_degrees.append({
        'node_id': node_id,
        'protein_id': protein_id,
        'degree': degree
    })

# Sort by degree (descending)
missing_node_degrees.sort(key=lambda x: x['degree'], reverse=True)

# Create DataFrame for analysis
missing_df = pd.DataFrame(missing_node_degrees)

print("\n=== MISSING NODE DEGREE ANALYSIS ===")
print(f"\nTotal missing nodes: {len(missing_nodes)}")
print(f"\nDegree statistics:")
print(missing_df['degree'].describe())

print(f"\nDegree distribution:")
print(missing_df['degree'].value_counts().sort_index())

print(f"\nTop 20 missing nodes by degree:")
print(missing_df.head(20).to_string(index=False))

print(f"\nNodes with degree 0 (isolated): {len(missing_df[missing_df['degree'] == 0])}")
print(f"Nodes with degree 1: {len(missing_df[missing_df['degree'] == 1])}")
print(f"Nodes with degree >= 5: {len(missing_df[missing_df['degree'] >= 5])}")
print(f"Nodes with degree >= 10: {len(missing_df[missing_df['degree'] >= 10])}")

# Save to CSV for further analysis
missing_df.to_csv('../processed_data/missing_nodes_degree_analysis.csv', index=False)
print(f"\nDetailed analysis saved to: ../processed_data/missing_nodes_degree_analysis.csv")
