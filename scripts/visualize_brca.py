import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

def visualize_ppi_network(edgelist_file, info_file, output_file='network_viz.png'):
    """
    Visualize a PPI network from edgelist and info files
    
    Args:
        edgelist_file: Path to edgelist (tab-separated: protein1, protein2, score)
        info_file: Path to protein info file
        output_file: Where to save the visualization
    """
    
    # Load edgelist
    print("Loading edgelist...")
    edgelist = pd.read_csv(edgelist_file, sep='\t', header=None,
                          names=['protein1', 'protein2', 'score'])
    
    # Load protein info
    print("Loading protein info...")
    info = pd.read_csv(info_file, sep='\t')
    protein_id_col = info.columns[0]  # First column is protein ID
    
    # Create mapping from protein ID to gene name
    id_to_name = dict(zip(info[protein_id_col], info['preferred_name']))
    
    # Build graph
    print("Building network...")
    G = nx.from_pandas_edgelist(
        edgelist,
        source='protein1',
        target='protein2',
        edge_attr='score'
    )
    
    print(f"\nNetwork Statistics:")
    print(f"  Nodes: {G.number_of_nodes()}")
    print(f"  Edges: {G.number_of_edges()}")
    print(f"  Average degree: {sum(dict(G.degree()).values()) / G.number_of_nodes():.2f}")
    
    # Get node labels (gene names)
    labels = {node: id_to_name.get(node, node.split('.')[-1][:8]) 
              for node in G.nodes()}
    
    # Calculate node sizes based on degree
    degrees = dict(G.degree())
    node_sizes = [degrees[node] * 50 for node in G.nodes()]
    
    # Find hub proteins
    top_hubs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:5]
    print(f"\nTop 5 Hub Proteins:")
    for node, degree in top_hubs:
        gene_name = id_to_name.get(node, 'Unknown')
        print(f"  {gene_name}: {degree} connections")
    
    # Create visualization
    print("\nCreating visualization...")
    plt.figure(figsize=(16, 12))
    
    # Use spring layout for better visualization
    pos = nx.spring_layout(G, k=0.5, iterations=50, seed=42)
    
    # Draw edges
    nx.draw_networkx_edges(
        G, pos,
        alpha=0.2,
        width=0.5,
        edge_color='gray'
    )
    
    # Draw nodes
    nx.draw_networkx_nodes(
        G, pos,
        node_size=node_sizes,
        node_color=list(degrees.values()),
        cmap='YlOrRd',
        alpha=0.8,
        edgecolors='black',
        linewidths=0.5
    )
    
    # Draw labels
    nx.draw_networkx_labels(
        G, pos,
        labels,
        font_size=8,
        font_weight='bold'
    )
    
    plt.title(f'Protein-Protein Interaction Network\n{G.number_of_nodes()} proteins, {G.number_of_edges()} interactions',
              fontsize=16, fontweight='bold')
    plt.axis('off')
    plt.tight_layout()
    
    # Save figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"\nVisualization saved to {output_file}")
    
    # Also create a degree distribution plot
    plt.figure(figsize=(10, 6))
    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)
    plt.hist(degree_sequence, bins=30, color='skyblue', edgecolor='black', alpha=0.7)
    plt.xlabel('Degree', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title('Degree Distribution', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    plt.savefig('degree_distribution.png', dpi=300, bbox_inches='tight')
    print("Degree distribution saved to degree_distribution.png")
    
    return G, id_to_name

# Usage
if __name__ == "__main__":
    G, id_to_name = visualize_ppi_network(
        edgelist_file='breast_cancer_ppi.txt',
        info_file='9606.protein.info.v12.0.txt',
        output_file='breast_cancer_network.png'
    )