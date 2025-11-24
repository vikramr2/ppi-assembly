#!/usr/bin/env python3
"""
Process BioGRID protein-protein interaction data to create:
1. Edge list CSV with node IDs
2. Metadata CSV mapping node IDs to gene names
"""

import csv
from pathlib import Path

def process_biogrid_data(input_file, output_dir):
    """
    Process BioGRID tab-delimited file and create edge list and metadata files.

    Args:
        input_file: Path to the BioGRID .tab3.txt file
        output_dir: Directory where output files will be saved
    """
    # Track unique genes and their information
    gene_map = {}  # biogrid_id -> gene_info
    edges = []  # list of (biogrid_id_a, biogrid_id_b) tuples

    print(f"Reading data from: {input_file}")

    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            # Extract BioGRID IDs (these will be our node IDs)
            biogrid_id_a = row['BioGRID ID Interactor A']
            biogrid_id_b = row['BioGRID ID Interactor B']

            # Extract gene names (Official Symbols)
            gene_name_a = row['Official Symbol Interactor A']
            gene_name_b = row['Official Symbol Interactor B']

            # Extract organism info for context
            organism_a = row['Organism Name Interactor A']
            organism_b = row['Organism Name Interactor B']

            # Extract Entrez Gene IDs for additional reference
            entrez_a = row['Entrez Gene Interactor A']
            entrez_b = row['Entrez Gene Interactor B']

            # Store gene information
            if biogrid_id_a not in gene_map:
                gene_map[biogrid_id_a] = {
                    'gene_name': gene_name_a,
                    'organism': organism_a,
                    'entrez_id': entrez_a
                }

            if biogrid_id_b not in gene_map:
                gene_map[biogrid_id_b] = {
                    'gene_name': gene_name_b,
                    'organism': organism_b,
                    'entrez_id': entrez_b
                }

            # Add edge (interaction)
            edges.append((biogrid_id_a, biogrid_id_b))

    print(f"Found {len(edges)} interactions")
    print(f"Found {len(gene_map)} unique genes")

    # Create output directory if it doesn't exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Write edge list CSV
    edge_file = output_dir / 'edgelist.csv'
    print(f"\nWriting edge list to: {edge_file}")
    with open(edge_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['source', 'target'])
        writer.writerows(edges)

    # Write metadata CSV
    metadata_file = output_dir / 'metadata.csv'
    print(f"Writing metadata to: {metadata_file}")
    with open(metadata_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['node_id', 'gene_name', 'entrez_id', 'organism'])

        for biogrid_id in sorted(gene_map.keys()):
            info = gene_map[biogrid_id]
            writer.writerow([
                biogrid_id,
                info['gene_name'],
                info['entrez_id'],
                info['organism']
            ])

    print("\nProcessing complete!")
    print(f"Edge list: {edge_file} ({len(edges)} edges)")
    print(f"Metadata: {metadata_file} ({len(gene_map)} nodes)")

if __name__ == '__main__':
    input_file = '/Users/vikram/Downloads/BIOGRID-GENE-4383915-5.0.251.DOWNLOADS/BIOGRID-GENE-4383915-5.0.251.tab3.txt'
    output_dir = '/Users/vikram/Documents/School/CS546/ppi-assembly'

    process_biogrid_data(input_file, output_dir)
