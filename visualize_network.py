#!/usr/bin/env python3
"""
Visualize the protein-protein interaction network with gene name labels.
Creates a beautiful, publication-quality network graph.
"""

import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import numpy as np

def load_network_data(edgelist_path, metadata_path):
    """
    Load the edgelist and metadata, creating a network with gene labels.

    Args:
        edgelist_path: Path to edgelist.csv
        metadata_path: Path to metadata.csv with node_id and gene_symbol

    Returns:
        NetworkX graph with gene symbols as labels
    """
    # Load edgelist
    edges = pd.read_csv(edgelist_path)
    print(f"Loaded {len(edges)} edges")

    # Load metadata
    metadata = pd.read_csv(metadata_path)
    print(f"Loaded metadata for {len(metadata)} genes")

    # Create mapping from node_id to gene_symbol
    node_to_gene = dict(zip(metadata['node_id'], metadata['gene_symbol']))

    # Create graph
    G = nx.from_pandas_edgelist(edges, source='source', target='target')

    # Remove self-loops for cleaner visualization
    G.remove_edges_from(nx.selfloop_edges(G))

    print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    return G, node_to_gene


def create_beautiful_visualization(G, node_to_gene, output_file='network_visualization.png'):
    """
    Create a beautiful network visualization with gene labels.

    Args:
        G: NetworkX graph
        node_to_gene: Dictionary mapping node IDs to gene symbols
        output_file: Output filename
    """
    # Set up the figure with high DPI for quality
    fig, ax = plt.subplots(figsize=(20, 20), dpi=150)
    fig.patch.set_facecolor('#f8f9fa')
    ax.set_facecolor('#ffffff')

    # Calculate node degrees for sizing
    degrees = dict(G.degree())
    max_degree = max(degrees.values()) if degrees else 1

    # Node sizes based on degree (connectivity)
    node_sizes = [300 + (degrees[node] / max_degree) * 1500 for node in G.nodes()]

    # Use spring layout for organic look (can try other layouts too)
    print("Computing network layout...")
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Create color map based on degree (hubs are darker)
    node_colors = [degrees[node] for node in G.nodes()]

    # Draw edges with transparency
    nx.draw_networkx_edges(
        G, pos,
        edge_color='#94a3b8',
        width=0.5,
        alpha=0.3,
        ax=ax
    )

    # Draw nodes with gradient colors
    nodes = nx.draw_networkx_nodes(
        G, pos,
        node_size=node_sizes,
        node_color=node_colors,
        cmap=plt.cm.viridis,
        alpha=0.9,
        edgecolors='#1e293b',
        linewidths=1.5,
        ax=ax
    )

    # Create gene symbol labels
    labels = {}
    for node in G.nodes():
        # Use gene symbol if available, otherwise use node ID
        labels[node] = node_to_gene.get(node, str(node))

    # Draw labels with nice formatting
    # For large networks, only label high-degree nodes
    if G.number_of_nodes() > 50:
        # Label only top 30% by degree
        threshold = sorted(degrees.values(), reverse=True)[int(len(degrees) * 0.3)]
        labels_to_draw = {k: v for k, v in labels.items() if degrees[k] >= threshold}
        print(f"Labeling {len(labels_to_draw)} high-degree nodes (degree >= {threshold})")
    else:
        labels_to_draw = labels

    nx.draw_networkx_labels(
        G, pos,
        labels=labels_to_draw,
        font_size=8,
        font_weight='bold',
        font_color='#1e293b',
        font_family='sans-serif',
        ax=ax
    )

    # Add colorbar to show degree scale
    sm = plt.cm.ScalarMappable(
        cmap=plt.cm.viridis,
        norm=plt.Normalize(vmin=min(node_colors), vmax=max(node_colors))
    )
    sm.set_array([])
    cbar = plt.colorbar(sm, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Node Degree (Number of Connections)', fontsize=12, fontweight='bold')

    # Add title and stats
    title = 'Protein-Protein Interaction Network'
    stats = f'{G.number_of_nodes()} proteins | {G.number_of_edges()} interactions'

    ax.text(
        0.5, 0.98, title,
        transform=ax.transAxes,
        fontsize=24,
        fontweight='bold',
        ha='center',
        va='top',
        color='#1e293b'
    )

    ax.text(
        0.5, 0.95, stats,
        transform=ax.transAxes,
        fontsize=14,
        ha='center',
        va='top',
        color='#475569',
        style='italic'
    )

    # Remove axes
    ax.axis('off')

    # Add a subtle border
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Save with tight layout
    plt.tight_layout()
    print(f"Saving visualization to {output_file}...")
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
    print(f"✓ Saved high-resolution visualization to {output_file}")

    plt.close()


def create_hub_analysis(G, node_to_gene, top_n=20):
    """
    Create a visualization highlighting network hubs (highly connected proteins).

    Args:
        G: NetworkX graph
        node_to_gene: Dictionary mapping node IDs to gene symbols
        top_n: Number of top hubs to highlight
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 10), dpi=150)
    fig.patch.set_facecolor('#f8f9fa')

    # Calculate degree centrality
    degrees = dict(G.degree())
    sorted_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)

    # --- Left plot: Network with hub highlighting ---
    ax1.set_facecolor('#ffffff')
    pos = nx.spring_layout(G, k=2, iterations=50, seed=42)

    # Identify hub nodes (top N)
    hub_nodes = [node for node, _ in sorted_nodes[:top_n]]
    hub_degrees = [degrees[node] for node in hub_nodes]

    # Draw regular nodes
    regular_nodes = [n for n in G.nodes() if n not in hub_nodes]
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=regular_nodes,
        node_size=200,
        node_color='#cbd5e1',
        alpha=0.6,
        edgecolors='#64748b',
        linewidths=0.5,
        ax=ax1
    )

    # Draw hub nodes (larger and colored)
    nx.draw_networkx_nodes(
        G, pos,
        nodelist=hub_nodes,
        node_size=[300 + d * 50 for d in hub_degrees],
        node_color=hub_degrees,
        cmap=plt.cm.Reds,
        alpha=0.9,
        edgecolors='#7f1d1d',
        linewidths=2,
        ax=ax1
    )

    # Draw edges
    nx.draw_networkx_edges(
        G, pos,
        edge_color='#94a3b8',
        width=0.3,
        alpha=0.2,
        ax=ax1
    )

    # Label only hub nodes
    hub_labels = {node: node_to_gene.get(node, str(node)) for node in hub_nodes}
    nx.draw_networkx_labels(
        G, pos,
        labels=hub_labels,
        font_size=9,
        font_weight='bold',
        font_color='#7f1d1d',
        ax=ax1
    )

    ax1.set_title('Network Hubs (Highly Connected Proteins)',
                  fontsize=16, fontweight='bold', pad=20, color='#1e293b')
    ax1.axis('off')

    # --- Right plot: Bar chart of top hubs ---
    ax2.set_facecolor('#ffffff')

    # Get top hubs with gene names
    top_hubs = [(node_to_gene.get(node, str(node)), degree)
                for node, degree in sorted_nodes[:top_n]]

    genes, degrees_list = zip(*top_hubs)
    y_pos = np.arange(len(genes))

    # Create horizontal bar chart
    bars = ax2.barh(y_pos, degrees_list, color=plt.cm.Reds(np.linspace(0.5, 0.9, len(genes))))

    # Customize bars
    for bar in bars:
        bar.set_edgecolor('#7f1d1d')
        bar.set_linewidth(1.5)

    ax2.set_yticks(y_pos)
    ax2.set_yticklabels(genes, fontsize=10, fontweight='bold')
    ax2.invert_yaxis()
    ax2.set_xlabel('Number of Interactions', fontsize=12, fontweight='bold')
    ax2.set_title(f'Top {top_n} Most Connected Proteins',
                  fontsize=16, fontweight='bold', pad=20, color='#1e293b')

    # Add value labels on bars
    for i, (bar, degree) in enumerate(zip(bars, degrees_list)):
        ax2.text(degree + max(degrees_list) * 0.01, i, str(degree),
                va='center', fontsize=9, fontweight='bold', color='#1e293b')

    # Style the plot
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(axis='x', alpha=0.3, linestyle='--')

    plt.tight_layout()
    output_file = 'network_hubs.png'
    print(f"Saving hub analysis to {output_file}...")
    plt.savefig(output_file, dpi=300, bbox_inches='tight', facecolor=fig.get_facecolor())
    print(f"✓ Saved hub analysis to {output_file}")

    plt.close()


def print_network_statistics(G, node_to_gene):
    """Print interesting network statistics."""
    print("\n" + "="*60)
    print("NETWORK STATISTICS")
    print("="*60)

    print(f"\nBasic Stats:")
    print(f"  Nodes (proteins): {G.number_of_nodes()}")
    print(f"  Edges (interactions): {G.number_of_edges()}")
    print(f"  Average degree: {sum(dict(G.degree()).values()) / G.number_of_nodes():.2f}")

    # Degree distribution
    degrees = dict(G.degree())
    print(f"\nDegree Distribution:")
    print(f"  Max degree: {max(degrees.values())}")
    print(f"  Min degree: {min(degrees.values())}")

    # Top hubs
    sorted_nodes = sorted(degrees.items(), key=lambda x: x[1], reverse=True)[:10]
    print(f"\nTop 10 Hub Proteins:")
    for i, (node, degree) in enumerate(sorted_nodes, 1):
        gene = node_to_gene.get(node, str(node))
        print(f"  {i:2d}. {gene:15s} ({degree} interactions)")

    # Connected components
    if not nx.is_connected(G):
        components = list(nx.connected_components(G))
        print(f"\nConnected Components: {len(components)}")
        print(f"  Largest component: {len(max(components, key=len))} nodes")
    else:
        print(f"\nNetwork is fully connected!")

    print("\n" + "="*60 + "\n")


if __name__ == '__main__':
    # File paths
    edgelist_path = 'processed_data/edgelist.csv'
    metadata_path = 'processed_data/metadata.csv'

    print("Loading network data...")
    G, node_to_gene = load_network_data(edgelist_path, metadata_path)

    # Print statistics
    print_network_statistics(G, node_to_gene)

    # Create main visualization
    print("Creating main network visualization...")
    create_beautiful_visualization(G, node_to_gene, 'network_visualization.png')

    # Create hub analysis
    print("\nCreating hub analysis...")
    create_hub_analysis(G, node_to_gene, top_n=20)

    print("\n✓ All visualizations complete!")
    print("  - network_visualization.png: Full network with gene labels")
    print("  - network_hubs.png: Hub analysis with top connected proteins")