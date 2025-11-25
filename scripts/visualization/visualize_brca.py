import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

def visualize_ppi_network(edgelist_file, info_file, output_file='network_viz.png',
                         layout='fa2', label_top_n=50):
    """
    Visualize a large PPI network (2000-3000 nodes) from edgelist and info files

    Args:
        edgelist_file: Path to edgelist CSV (columns: source, target)
        info_file: Path to protein info CSV file
        output_file: Where to save the visualization
        layout: Layout algorithm ('fa2' for ForceAtlas2, 'spring', or 'kamada_kawai')
        label_top_n: Number of top hub proteins to label (default 50)
    """

    # Load edgelist
    print("Loading edgelist...")
    edgelist = pd.read_csv(edgelist_file)

    # Load protein info
    print("Loading protein info...")
    info = pd.read_csv(info_file)

    # Create mapping from node_id to gene name
    id_to_name = dict(zip(info['node_id'], info['preferred_name']))

    # Build graph
    print("Building network...")
    G = nx.from_pandas_edgelist(
        edgelist,
        source='source',
        target='target'
    )

    print(f"\nNetwork Statistics:")
    print(f"  Nodes: {G.number_of_nodes()}")
    print(f"  Edges: {G.number_of_edges()}")
    print(f"  Average degree: {sum(dict(G.degree()).values()) / G.number_of_nodes():.2f}")
    print(f"  Density: {nx.density(G):.6f}")

    # Calculate node sizes based on degree
    degrees = dict(G.degree())

    # Find hub proteins
    top_hubs = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:20]
    print(f"\nTop 20 Hub Proteins:")
    for node, degree in top_hubs:
        gene_name = id_to_name.get(node, 'Unknown')
        print(f"  {gene_name}: {degree} connections")

    # Create visualization
    print(f"\nCreating visualization using {layout} layout...")
    print("This may take a few minutes for large graphs...")

    plt.figure(figsize=(24, 24))

    # Choose layout algorithm based on size
    pos = None
    if layout == 'fa2':
        try:
            from fa2 import ForceAtlas2
            forceatlas2 = ForceAtlas2(
                outboundAttractionDistribution=True,
                linLogMode=False,
                adjustSizes=False,
                edgeWeightInfluence=1.0,
                jitterTolerance=1.0,
                barnesHutOptimize=True,
                barnesHutTheta=1.2,
                multiThreaded=False,
                scalingRatio=2.0,
                strongGravityMode=False,
                gravity=1.0,
                verbose=True
            )
            print("Computing ForceAtlas2 layout...")
            pos = forceatlas2.forceatlas2_networkx_layout(G, pos=None, iterations=2000)
        except ImportError:
            print("ForceAtlas2 not available, falling back to spring layout")
            layout = 'spring'

    if layout == 'spring' or pos is None:
        print("Computing spring layout with sparse optimization...")
        pos = nx.spring_layout(G, k=2/np.sqrt(G.number_of_nodes()),
                              iterations=50, seed=42)
    elif layout == 'kamada_kawai':
        print("Computing Kamada-Kawai layout...")
        pos = nx.kamada_kawai_layout(G)

    # Normalize node sizes for large networks
    max_degree = max(degrees.values())
    node_sizes = [20 + (degrees[node] / max_degree) * 200 for node in G.nodes()]

    # Color nodes by degree with better contrast
    node_colors = [degrees[node] for node in G.nodes()]

    # Draw edges with very low opacity for dense networks
    print("Drawing edges...")
    nx.draw_networkx_edges(
        G, pos,
        alpha=0.25,
        width=0.2,
        edge_color='gray'
    )

    # Draw nodes
    print("Drawing nodes...")
    nx.draw_networkx_nodes(
        G, pos,
        node_size=node_sizes,
        node_color=node_colors,
        cmap='viridis',
        alpha=0.7,
        edgecolors='white',
        linewidths=0.3
    )

    # Only label top hub proteins to avoid clutter
    print(f"Adding labels for top {label_top_n} hubs...")
    top_nodes = [node for node, _ in sorted(degrees.items(), key=lambda x: x[1],
                                            reverse=True)[:label_top_n]]
    labels = {node: id_to_name.get(node, str(node)) for node in top_nodes}

    nx.draw_networkx_labels(
        G, pos,
        labels,
        font_size=6,
        font_weight='bold',
        font_color='black',
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white',
                 edgecolor='none', alpha=0.7)
    )
    
    plt.title(f'Breast Cancer PPI Network\n{G.number_of_nodes()} proteins, {G.number_of_edges()} interactions',
              fontsize=20, fontweight='bold', pad=20)
    plt.axis('off')
    plt.tight_layout()

    # Save figure
    print("Saving visualization...")
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"Visualization saved to {output_file}")
    plt.close()

    # Create degree distribution plot
    print("Creating degree distribution plot...")
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    degree_sequence = sorted([d for n, d in G.degree()], reverse=True)

    # Linear scale
    ax1.hist(degree_sequence, bins=50, color='skyblue', edgecolor='black', alpha=0.7)
    ax1.set_xlabel('Degree', fontsize=12)
    ax1.set_ylabel('Frequency', fontsize=12)
    ax1.set_title('Degree Distribution (Linear)', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)

    # Log-log scale
    degree_counts = {}
    for degree in degree_sequence:
        degree_counts[degree] = degree_counts.get(degree, 0) + 1

    degrees_x = sorted(degree_counts.keys())
    counts_y = [degree_counts[d] for d in degrees_x]

    ax2.loglog(degrees_x, counts_y, 'o', markersize=8, alpha=0.6, color='coral')
    ax2.set_xlabel('Degree (log scale)', fontsize=12)
    ax2.set_ylabel('Frequency (log scale)', fontsize=12)
    ax2.set_title('Degree Distribution (Log-Log)', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, which='both')

    plt.tight_layout()
    plt.savefig('degree_distribution.png', dpi=300, bbox_inches='tight')
    print("Degree distribution saved to degree_distribution.png")
    plt.close()

    return G, id_to_name

# Usage
if __name__ == "__main__":
    import sys

    # Default files
    edgelist_file = '../processed_data/brca_ppi_edgelist_cleaned.csv'
    info_file = '../processed_data/brca_protein_info.csv'
    output_file = 'brca_network_visualization.png'

    # Parse command line arguments if provided
    if len(sys.argv) > 1:
        edgelist_file = sys.argv[1]
    if len(sys.argv) > 2:
        info_file = sys.argv[2]
    if len(sys.argv) > 3:
        output_file = sys.argv[3]

    print("=" * 60)
    print("BRCA PPI Network Visualization")
    print("=" * 60)

    G, id_to_name = visualize_ppi_network(
        edgelist_file=edgelist_file,
        info_file=info_file,
        output_file=output_file,
        layout='spring',  # Use 'fa2' if you have ForceAtlas2 installed
        label_top_n=50
    )

    print("\n" + "=" * 60)
    print("Visualization complete!")
    print("=" * 60)