# Protein-Protein Interaction (PPI) Network Analysis

## Overview

This dataset contains a breast cancer-focused protein-protein interaction (PPI) network derived from the STRING database. The network was constructed by starting with a seed set of breast cancer-related proteins and expanding one step via breadth-first search (BFS) to include direct interaction partners.

## Data Collection Pipeline

1. **Source**: Downloaded PPI data from [STRING database](https://string-db.org/)
2. **Seed Set**: Started with breast cancer-associated proteins
3. **Network Expansion**: Performed 1-step BFS to include direct neighbors
4. **Sequence Retrieval**: Fetched amino acid sequences using Ensembl REST API

## Processed Data Files

### `processed_data/brca_protein_info.csv`
- **Size**: 1.1 MB
- **Entries**: 2,394 proteins
- **Columns**:
  - `protein_id`: Ensembl protein ID (format: 9606.ENSP########)
  - `preferred_name`: Gene symbol (e.g., ARF5, BRCA1)
  - `protein_size`: Number of amino acids
  - `annotation`: Functional description from STRING
  - `node_id`: Integer node identifier (0-2393)

### `processed_data/brca_ppi_edgelist_cleaned.csv`
- **Size**: 526 KB
- **Entries**: 53,363 unique edges
- **Format**: Undirected edge list (source ≤ target)
- **Columns**:
  - `source`: Node ID of first protein
  - `target`: Node ID of second protein
- **Notes**: Duplicate and reverse edges removed to represent undirected graph

### `processed_data/node_to_sequence.json`
- **Size**: 1.5 MB
- **Format**: JSON mapping
- **Structure**: `{node_id: amino_acid_sequence}`
- **Coverage**: Amino acid sequences for proteins in the network
- **Source**: Retrieved from Ensembl REST API

## Network Statistics

- **Proteins**: 2,394
- **Interactions**: 53,363 (undirected)
- **Organism**: Human (Homo sapiens, NCBI taxonomy ID 9606)

## Usage

The node IDs (0-2393) are used consistently across all files:
- Use `brca_protein_info.csv` to map node IDs to protein identifiers and metadata
- Use `brca_ppi_edgelist_cleaned.csv` for network topology
- Use `node_to_sequence.json` to access amino acid sequences for each protein

## Scripts

- `fetch_protein.py`: Retrieves protein sequences from Ensembl API
- `check_parallel_reverse_edges.py`: Cleans edge list by removing duplicates and reverse edges
