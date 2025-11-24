#!/usr/bin/env python3
"""
Create an interactive HTML visualization of the PPI network using Plotly.
This allows zooming, panning, and hovering to see gene information.
"""

import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from plotly.offline import plot
import numpy as np

def load_network_data(edgelist_path, metadata_path):
    """Load network and metadata."""
    edges = pd.read_csv(edgelist_path)
    metadata = pd.read_csv(metadata_path)

    # Create mapping from node_id to gene info
    node_info = {}
    for _, row in metadata.iterrows():
        summary = row.get('summary', 'No summary available')
        if pd.notna(summary) and isinstance(summary, str):
            summary = summary[:200] + '...' if len(summary) > 200 else summary
        else:
            summary = 'No summary available'

        node_info[row['node_id']] = {
            'gene_symbol': row['gene_symbol'],
            'gene_name': row['gene_name'],
            'chromosome': row.get('chromosome', 'N/A'),
            'summary': summary
        }

    # Create graph
    G = nx.from_pandas_edgelist(edges, source='source', target='target')
    G.remove_edges_from(nx.selfloop_edges(G))

    return G, node_info


def create_interactive_plot(G, node_info, output_file='network_interactive.html'):
    """Create interactive Plotly visualization."""

    print("Computing layout for interactive visualization...")

    # For star networks, use a custom radial layout
    degrees = dict(G.degree())
    hub_node = max(degrees.items(), key=lambda x: x[1])[0]

    if degrees[hub_node] > G.number_of_nodes() * 0.8:
        # Star network - use radial layout
        print(f"Detected star network with hub: {node_info.get(hub_node, {}).get('gene_symbol', hub_node)}")
        pos = create_radial_layout(G, hub_node)
    else:
        # Regular network - use spring layout
        pos = nx.spring_layout(G, k=1, iterations=50, seed=42)

    # Extract edge coordinates
    edge_x = []
    edge_y = []
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])

    edge_trace = go.Scatter(
        x=edge_x, y=edge_y,
        line=dict(width=0.5, color='#888'),
        hoverinfo='none',
        mode='lines',
        opacity=0.3
    )

    # Extract node coordinates and info
    node_x = []
    node_y = []
    node_text = []
    node_colors = []
    node_sizes = []

    max_degree = max(degrees.values())

    for node in G.nodes():
        x, y = pos[node]
        node_x.append(x)
        node_y.append(y)

        # Get gene info
        info = node_info.get(node, {})
        gene_symbol = info.get('gene_symbol', str(node))
        gene_name = info.get('gene_name', 'Unknown')
        chromosome = info.get('chromosome', 'N/A')
        summary = info.get('summary', 'No summary')
        degree = degrees[node]

        # Create hover text
        hover_text = f"<b>{gene_symbol}</b><br>"
        hover_text += f"Full name: {gene_name}<br>"
        hover_text += f"Chromosome: {chromosome}<br>"
        hover_text += f"Interactions: {degree}<br>"
        hover_text += f"<br>{summary}"

        node_text.append(hover_text)
        node_colors.append(degree)

        # Size based on degree
        size = 10 + (degree / max_degree) * 40
        node_sizes.append(size)

    node_trace = go.Scatter(
        x=node_x, y=node_y,
        mode='markers+text',
        hoverinfo='text',
        text=[node_info.get(node, {}).get('gene_symbol', str(node)) for node in G.nodes()],
        textposition="top center",
        textfont=dict(size=8, color='#1e293b'),
        hovertext=node_text,
        marker=dict(
            showscale=True,
            colorscale='Viridis',
            size=node_sizes,
            color=node_colors,
            colorbar=dict(
                thickness=15,
                title='Node Degree',
                xanchor='left',
                titleside='right'
            ),
            line=dict(width=2, color='#1e293b')
        )
    )

    # Create figure
    fig = go.Figure(data=[edge_trace, node_trace],
                   layout=go.Layout(
                       title=dict(
                           text='<b>Interactive Protein-Protein Interaction Network</b><br>' +
                                f'<i>{G.number_of_nodes()} proteins | {G.number_of_edges()} interactions</i>',
                           x=0.5,
                           xanchor='center',
                           font=dict(size=20)
                       ),
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20, l=5, r=5, t=100),
                       annotations=[dict(
                           text="Hover over nodes to see protein details | Scroll to zoom | Drag to pan",
                           showarrow=False,
                           xref="paper", yref="paper",
                           x=0.5, y=-0.05,
                           xanchor='center', yanchor='top',
                           font=dict(size=12, color='#64748b')
                       )],
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False),
                       plot_bgcolor='#f8f9fa',
                       paper_bgcolor='#ffffff'
                   ))

    print(f"Saving interactive visualization to {output_file}...")
    plot(fig, filename=output_file, auto_open=False)
    print(f"✓ Saved interactive visualization to {output_file}")
    print(f"  Open {output_file} in your browser to explore the network!")


def create_radial_layout(G, hub_node):
    """Create a radial layout for star networks with hub at center."""
    pos = {}

    # Place hub at center
    pos[hub_node] = (0, 0)

    # Place other nodes in concentric circles
    other_nodes = [n for n in G.nodes() if n != hub_node]
    n_nodes = len(other_nodes)

    # Calculate number of rings needed
    nodes_per_ring = 30
    n_rings = max(1, (n_nodes + nodes_per_ring - 1) // nodes_per_ring)

    node_idx = 0
    for ring in range(n_rings):
        radius = (ring + 1) * 1.5
        nodes_in_ring = min(nodes_per_ring, n_nodes - node_idx)

        for i in range(nodes_in_ring):
            angle = 2 * np.pi * i / nodes_in_ring
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            pos[other_nodes[node_idx]] = (x, y)
            node_idx += 1

    return pos


if __name__ == '__main__':
    # File paths
    edgelist_path = 'processed_data/edgelist.csv'
    metadata_path = 'processed_data/metadata.csv'

    print("Loading network data...")
    G, node_info = load_network_data(edgelist_path, metadata_path)

    print(f"Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # Create interactive visualization
    create_interactive_plot(G, node_info)

    print("\n✓ Interactive visualization complete!")