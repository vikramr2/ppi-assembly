import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load the Leiden clustering
print("Loading Leiden clustering...")
leiden_df = pd.read_csv('../../processed_data/clustering/disjoint/brca_ppi_leiden.csv')
print(f"Loaded {len(leiden_df)} node assignments")

# Calculate cluster sizes
cluster_sizes = leiden_df['cluster_id'].value_counts().sort_index()
print(f"\nFound {len(cluster_sizes)} clusters")

print(f"\nCluster size statistics:")
print(cluster_sizes.describe())

# Create histogram
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Histogram 1: All clusters
ax1.hist(cluster_sizes.values, bins=50, edgecolor='black', alpha=0.7, color='steelblue')
ax1.set_xlabel('Cluster Size (number of proteins)', fontsize=12)
ax1.set_ylabel('Number of Clusters', fontsize=12)
ax1.set_title(f'Leiden Cluster Size Distribution\n({len(cluster_sizes)} clusters total)', fontsize=13)
ax1.grid(axis='y', alpha=0.3)

# Add statistics text
stats_text = f'Mean: {cluster_sizes.mean():.1f}\nMedian: {cluster_sizes.median():.1f}\nMax: {cluster_sizes.max()}'
ax1.text(0.65, 0.95, stats_text, transform=ax1.transAxes,
         verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Histogram 2: Log scale for better visualization
ax2.hist(cluster_sizes.values, bins=50, edgecolor='black', alpha=0.7, color='coral')
ax2.set_xlabel('Cluster Size (number of proteins)', fontsize=12)
ax2.set_ylabel('Number of Clusters (log scale)', fontsize=12)
ax2.set_title('Leiden Cluster Size Distribution (Log Scale)', fontsize=13)
ax2.set_yscale('log')
ax2.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.savefig('../../processed_data/clustering/disjoint/leiden_cluster_histogram.png', dpi=300, bbox_inches='tight')
print(f"\nHistogram saved to: ../processed_data/clustering/disjoint/leiden_cluster_histogram.png")

# Print detailed cluster information
print(f"\n=== CLUSTER SIZE BREAKDOWN ===")
print(f"Clusters with 1 protein: {len(cluster_sizes[cluster_sizes == 1])}")
print(f"Clusters with 2-10 proteins: {len(cluster_sizes[(cluster_sizes >= 2) & (cluster_sizes <= 10)])}")
print(f"Clusters with 11-50 proteins: {len(cluster_sizes[(cluster_sizes >= 11) & (cluster_sizes <= 50)])}")
print(f"Clusters with 51-100 proteins: {len(cluster_sizes[(cluster_sizes >= 51) & (cluster_sizes <= 100)])}")
print(f"Clusters with >100 proteins: {len(cluster_sizes[cluster_sizes > 100])}")

print(f"\nTop 10 largest clusters:")
for idx, (cluster_id, size) in enumerate(cluster_sizes.nlargest(10).items(), 1):
    print(f"  {idx}. Cluster {cluster_id}: {size} proteins")

# Save cluster size summary
cluster_summary = pd.DataFrame({
    'cluster_id': cluster_sizes.index,
    'cluster_size': cluster_sizes.values
}).sort_values('cluster_size', ascending=False)
cluster_summary.to_csv('../../processed_data/clustering/disjoint/leiden_cluster_sizes.csv', index=False)
print(f"\nCluster sizes saved to: ../processed_data/clustering/disjoint/leiden_cluster_sizes.csv")

plt.show()
