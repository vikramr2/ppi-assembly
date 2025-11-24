#!/usr/bin/env python3
"""
Fetch gene properties from NCBI using Entrez IDs.
Demonstrates what information can be retrieved for each gene.
"""

from Bio import Entrez
import csv
import json
import time
from pathlib import Path

# IMPORTANT: Always provide your email when using NCBI E-utilities
Entrez.email = "your.email@example.com"  # Change this to your email

def fetch_gene_summary(entrez_id):
    """
    Fetch gene summary information using Entrez ID.

    Args:
        entrez_id: NCBI Entrez Gene ID (as string or int)

    Returns:
        Dictionary with gene summary information
    """
    try:
        # Use esummary to get summary information
        handle = Entrez.esummary(db="gene", id=str(entrez_id))
        record = Entrez.read(handle)
        handle.close()

        if 'DocumentSummarySet' in record and 'DocumentSummary' in record['DocumentSummarySet']:
            summary = record['DocumentSummarySet']['DocumentSummary'][0]
            return summary
        return None
    except Exception as e:
        print(f"Error fetching {entrez_id}: {e}")
        return None

def fetch_gene_detailed(entrez_id):
    """
    Fetch detailed gene information in XML format.

    Args:
        entrez_id: NCBI Entrez Gene ID

    Returns:
        Raw XML data with comprehensive gene information
    """
    try:
        handle = Entrez.efetch(db="gene", id=str(entrez_id), retmode="xml")
        data = Entrez.read(handle)
        handle.close()
        return data
    except Exception as e:
        print(f"Error fetching detailed info for {entrez_id}: {e}")
        return None

def extract_properties(summary, detailed=None):
    """
    Extract useful properties from gene summary.

    Args:
        summary: Gene summary from esummary
        detailed: Optional detailed data from efetch

    Returns:
        Dictionary of extracted properties
    """
    if not summary:
        return {}

    properties = {
        'entrez_id': summary.get('Id', ''),
        'gene_symbol': summary.get('NomenclatureSymbol', summary.get('Name', '')),
        'gene_name': summary.get('Description', ''),
        'organism': summary.get('Organism', {}).get('ScientificName', ''),
        'chromosome': summary.get('Chromosome', ''),
        'map_location': summary.get('MapLocation', ''),
        'gene_type': summary.get('GeneType', ''),
        'summary': summary.get('Summary', ''),
        # These are already strings from NCBI, not lists
        'other_aliases': summary.get('OtherAliases', ''),
        'other_designations': summary.get('OtherDesignations', ''),
    }

    # Extract genomic info if available
    if 'GenomicInfo' in summary and summary['GenomicInfo']:
        genomic = summary['GenomicInfo'][0]
        properties['chr_accession'] = genomic.get('ChrAccVer', '')
        properties['chr_start'] = genomic.get('ChrStart', '')
        properties['chr_stop'] = genomic.get('ChrStop', '')

    return properties

def process_metadata_file(metadata_file, output_file, sample_size=245):
    """
    Read metadata CSV and fetch properties for genes.

    Args:
        metadata_file: Path to metadata.csv
        output_file: Path to save enriched data
        sample_size: Number of genes to process (use None for all)
    """
    print(f"Reading metadata from: {metadata_file}")

    genes = []
    with open(metadata_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row['entrez_id'] and row['entrez_id'] != '-':
                genes.append(row)

    print(f"Found {len(genes)} genes with Entrez IDs")

    # Process sample
    if sample_size:
        genes = genes[:sample_size]
        print(f"Processing first {sample_size} genes as example...")

    enriched_data = []

    for i, gene in enumerate(genes, 1):
        entrez_id = gene['entrez_id']
        print(f"\n[{i}/{len(genes)}] Fetching data for {gene['gene_name']} (Entrez: {entrez_id})...")

        # Fetch summary
        summary = fetch_gene_summary(entrez_id)
        if summary:
            properties = extract_properties(summary)
            enriched_data.append(properties)

            # Print sample properties
            print(f"  Symbol: {properties.get('gene_symbol', 'N/A')}")
            print(f"  Type: {properties.get('gene_type', 'N/A')}")
            print(f"  Location: Chr {properties.get('chromosome', 'N/A')} ({properties.get('map_location', 'N/A')})")
            if properties.get('summary'):
                print(f"  Summary: {properties['summary'][:100]}...")

        # Respect NCBI rate limit (3 requests/second without API key)
        time.sleep(0.35)

    # Save enriched data
    if enriched_data:
        output_path = Path(output_file)

        # Save as CSV
        csv_file = output_path.with_suffix('.csv')
        print(f"\nSaving enriched data to: {csv_file}")
        with open(csv_file, 'w', newline='', encoding='utf-8') as f:
            fieldnames = enriched_data[0].keys()
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(enriched_data)

        # Save as JSON for better structure
        json_file = output_path.with_suffix('.json')
        print(f"Saving enriched data to: {json_file}")
        with open(json_file, 'w', encoding='utf-8') as f:
            json.dump(enriched_data, f, indent=2)

        print(f"\nProcessed {len(enriched_data)} genes successfully!")
        return enriched_data

    return []

def demonstrate_detailed_fetch(entrez_id="59272"):
    """
    Demonstrate fetching detailed information for a single gene.
    Default: ACE2 (Entrez ID: 59272)
    """
    print(f"\n{'='*60}")
    print(f"DETAILED INFORMATION EXAMPLE - Entrez ID: {entrez_id}")
    print(f"{'='*60}\n")

    summary = fetch_gene_summary(entrez_id)
    if summary:
        print("Available fields in summary:")
        for key in sorted(summary.keys()):
            print(f"  - {key}")

        print("\n" + "="*60)
        properties = extract_properties(summary)
        print("\nExtracted Properties:")
        for key, value in properties.items():
            if value:
                print(f"  {key}: {value}")

if __name__ == '__main__':
    # Example 1: Show what detailed information looks like for ACE2
    demonstrate_detailed_fetch("59272")  # ACE2

    # Example 2: Process metadata file
    metadata_file = '/Users/vikram/Documents/School/CS546/ppi-assembly/metadata.csv'
    output_file = '/Users/vikram/Documents/School/CS546/ppi-assembly/enriched_metadata'

    print("\n" + "="*60)
    print("PROCESSING METADATA FILE")
    print("="*60)

    # Process only first 10 genes as example (change sample_size=None to process all)
    enriched_data = process_metadata_file(metadata_file, output_file, sample_size=245)
