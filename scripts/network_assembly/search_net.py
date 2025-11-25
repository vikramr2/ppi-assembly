import pandas as pd
import networkx as nx

def expand_network_bfs(seed_genes, ppi_df, info_df, min_score=700, max_hops=1, max_nodes=500):
    """
    Expand network from seed genes using BFS
    
    Args:
        seed_genes: List of gene names to start from
        ppi_df: DataFrame with all PPIs
        info_df: DataFrame with protein info
        min_score: Minimum interaction score
        max_hops: How many steps away from seed genes (1 = direct neighbors, 2 = neighbors of neighbors)
        max_nodes: Maximum number of nodes to include (to keep network manageable)
    
    Returns:
        DataFrame with expanded edgelist
    """
    
    protein_id_column = info_df.columns[0]
    
    # Get seed protein IDs
    seed_proteins = info_df[info_df['preferred_name'].isin(seed_genes)]
    seed_ids = set(seed_proteins[protein_id_column].values)
    
    print(f"Starting with {len(seed_ids)} seed proteins")
    
    # Filter high-confidence interactions
    high_conf_ppi = ppi_df[ppi_df['combined_score'] > min_score]
    
    # Build full graph for BFS
    G = nx.from_pandas_edgelist(
        high_conf_ppi,
        source='protein1',
        target='protein2',
        edge_attr='combined_score'
    )
    
    print(f"Full high-confidence network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Perform BFS to find neighbors
    expanded_ids = set(seed_ids)
    
    for hop in range(max_hops):
        print(f"\nHop {hop + 1}...")
        new_neighbors = set()
        
        for protein in expanded_ids:
            if protein in G:
                neighbors = set(G.neighbors(protein))
                new_neighbors.update(neighbors)
        
        expanded_ids.update(new_neighbors)
        print(f"  Total nodes after hop {hop + 1}: {len(expanded_ids)}")
        
        if len(expanded_ids) > max_nodes:
            print(f"  Reached max_nodes limit ({max_nodes})")
            break
    
    # If we exceeded max_nodes, keep only highest degree nodes
    if len(expanded_ids) > max_nodes:
        print(f"\nTrimming from {len(expanded_ids)} to {max_nodes} nodes...")
        
        # Calculate degree for all expanded nodes
        subgraph = G.subgraph(expanded_ids)
        degrees = dict(subgraph.degree())
        
        # Keep seed genes + highest degree nodes
        top_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)
        expanded_ids = seed_ids.copy()
        
        for node, degree in top_nodes:
            if len(expanded_ids) >= max_nodes:
                break
            expanded_ids.add(node)
    
    print(f"\nFinal network size: {len(expanded_ids)} nodes")
    
    # Extract subnetwork
    expanded_ppi = high_conf_ppi[
        (high_conf_ppi['protein1'].isin(expanded_ids)) & 
        (high_conf_ppi['protein2'].isin(expanded_ids))
    ]
    
    print(f"Final network: {len(expanded_ppi)} interactions")
    
    return expanded_ppi, expanded_ids


# ============ Main Script ============

# Load data
print("Loading data...")
ppi = pd.read_csv("9606.protein.links.v12.0.txt", sep=' ')
info = pd.read_csv("9606.protein.info.v12.0.txt", sep='\t')

print(f"Total interactions: {len(ppi)}")

# Seed genes
breast_cancer_genes = [
    'BRCA1', 'BRCA2', 'TP53', 'PTEN', 'CDH1', 'STK11',
    'ATM', 'CHEK2', 'PALB2', 'BARD1', 'BRIP1', 'RAD51C', 'RAD51D',
    'ESR1', 'ESR2', 'PGR', 'AR',
    'ERBB2', 'EGFR', 'ERBB3',
    'PIK3CA', 'PIK3R1', 'AKT1', 'AKT2', 'MTOR',
    'CCND1', 'CDK4', 'CDK6', 'CDKN2A', 'RB1',
    'GATA3', 'FOXA1', 'MYC', 'MAP3K1', 'ARID1A'
]

# Expand network with BFS
# Try different parameters:
# - max_hops=1: Direct neighbors only (medium network)
# - max_hops=2: Neighbors of neighbors (large network)
# - max_nodes: Controls final size

expanded_network, expanded_ids = expand_network_bfs(
    seed_genes=breast_cancer_genes,
    ppi_df=ppi,
    info_df=info,
    min_score=700,
    max_hops=1,        # Change to 2 for even larger network
    max_nodes=5000      # Adjust this for network size
)

# Save expanded network
expanded_network.to_csv("breast_cancer_ppi.txt", 
                       columns=['protein1', 'protein2', 'combined_score'],
                       sep='\t', index=False, header=False)

print("\nSaved to breast_cancer_ppi.txt")

# Print some stats about seed vs. expanded
protein_id_column = info.columns[0]
id_to_name = dict(zip(info[protein_id_column], info['preferred_name']))

print("\n=== Network Composition ===")
seed_ids = set(info[info['preferred_name'].isin(breast_cancer_genes)][protein_id_column].values)
new_proteins = expanded_ids - seed_ids

print(f"Seed proteins: {len(seed_ids)}")
print(f"New proteins added: {len(new_proteins)}")
print(f"\nExample new proteins:")
for protein_id in list(new_proteins)[:10]:
    gene_name = id_to_name.get(protein_id, 'Unknown')
    print(f"  {gene_name}")
