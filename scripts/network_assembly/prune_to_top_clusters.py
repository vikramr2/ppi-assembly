import pandas as pd

# Load the Leiden clustering
print("Loading Leiden clustering...")
leiden_df = pd.read_csv('../../processed_data/clustering/disjoint/brca_ppi_leiden.csv')
print(f"Loaded {len(leiden_df)} node assignments")

# Calculate cluster sizes
cluster_sizes = leiden_df['cluster_id'].value_counts()
print(f"\nFound {len(cluster_sizes)} clusters")
print(f"Total nodes: {len(leiden_df)}")

# Get top 5 largest clusters
top_5_clusters = cluster_sizes.nlargest(5)
print(f"\n=== TOP 5 LARGEST CLUSTERS ===")
for idx, (cluster_id, size) in enumerate(top_5_clusters.items(), 1):
    print(f"{idx}. Cluster {cluster_id}: {size} proteins ({size/len(leiden_df)*100:.1f}%)")

top_5_cluster_ids = set(top_5_clusters.index)

# Filter clustering to only top 5 clusters
leiden_pruned = leiden_df[leiden_df['cluster_id'].isin(top_5_cluster_ids)].copy()
print(f"\nPruned clustering: {len(leiden_pruned)} nodes ({len(leiden_pruned)/len(leiden_df)*100:.1f}% of original)")

# Get the set of node IDs to keep
nodes_to_keep = set(leiden_pruned['node_id'])

# Load and prune the edge list
print("\nLoading edge list...")
edgelist_df = pd.read_csv('../../processed_data/network/brca_ppi_edgelist_cleaned.csv')
print(f"Original edge list: {len(edgelist_df)} edges")

# Keep only edges where both nodes are in top 5 clusters
edgelist_pruned = edgelist_df[
    edgelist_df['source'].isin(nodes_to_keep) &
    edgelist_df['target'].isin(nodes_to_keep)
].copy()
print(f"Pruned edge list: {len(edgelist_pruned)} edges ({len(edgelist_pruned)/len(edgelist_df)*100:.1f}% of original)")

# Load and prune protein info
print("\nLoading protein info...")
protein_info_df = pd.read_csv('../../processed_data/metadata/brca_protein_info.csv')
print(f"Original protein info: {len(protein_info_df)} proteins")

protein_info_pruned = protein_info_df[protein_info_df['node_id'].isin(nodes_to_keep)].copy()
print(f"Pruned protein info: {len(protein_info_pruned)} proteins")

# Save pruned files
output_dir = '../../processed_data/clustering/disjoint/'

print("\n=== SAVING PRUNED FILES ===")

# Save pruned clustering
leiden_pruned.to_csv(f'{output_dir}brca_ppi_leiden_top5.csv', index=False)
print(f"Saved: {output_dir}brca_ppi_leiden_top5.csv")

# Save pruned edge list
edgelist_pruned.to_csv(f'{output_dir}brca_ppi_edgelist_top5.csv', index=False)
print(f"Saved: {output_dir}brca_ppi_edgelist_top5.csv")

# Save pruned protein info
protein_info_pruned.to_csv(f'{output_dir}brca_protein_info_top5.csv', index=False)
print(f"Saved: {output_dir}brca_protein_info_top5.csv")

# Create a summary file
summary = {
    'metric': [
        'Total clusters (original)',
        'Clusters kept',
        'Total nodes (original)',
        'Nodes in top 5 clusters',
        'Percentage of nodes kept',
        'Total edges (original)',
        'Edges in top 5 clusters',
        'Percentage of edges kept',
        'Cluster 1 ID',
        'Cluster 1 size',
        'Cluster 2 ID',
        'Cluster 2 size',
        'Cluster 3 ID',
        'Cluster 3 size',
        'Cluster 4 ID',
        'Cluster 4 size',
        'Cluster 5 ID',
        'Cluster 5 size',
    ],
    'value': [
        len(cluster_sizes),
        5,
        len(leiden_df),
        len(leiden_pruned),
        f"{len(leiden_pruned)/len(leiden_df)*100:.2f}%",
        len(edgelist_df),
        len(edgelist_pruned),
        f"{len(edgelist_pruned)/len(edgelist_df)*100:.2f}%",
        top_5_clusters.index[0],
        top_5_clusters.values[0],
        top_5_clusters.index[1],
        top_5_clusters.values[1],
        top_5_clusters.index[2],
        top_5_clusters.values[2],
        top_5_clusters.index[3],
        top_5_clusters.values[3],
        top_5_clusters.index[4],
        top_5_clusters.values[4],
    ]
}

summary_df = pd.DataFrame(summary)
summary_df.to_csv(f'{output_dir}top5_pruning_summary.csv', index=False)
print(f"Saved: {output_dir}top5_pruning_summary.csv")

print("\n=== SUMMARY ===")
print(f"Original network: {len(leiden_df)} nodes, {len(edgelist_df)} edges, {len(cluster_sizes)} clusters")
print(f"Pruned network: {len(leiden_pruned)} nodes, {len(edgelist_pruned)} edges, 5 clusters")
print(f"Reduction: {len(leiden_df) - len(leiden_pruned)} nodes removed, {len(edgelist_df) - len(edgelist_pruned)} edges removed")
