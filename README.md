# Protein-Protein Interaction (PPI) Network Analysis

This project processes BioGRID protein-protein interaction data and enriches it with gene properties from NCBI Entrez.

## Files Generated

### 1. Core Network Files
- **`edgelist.csv`** - Network edge list (354 interactions)
  - Columns: `source`, `target` (both are BioGRID node IDs)

- **`metadata.csv`** - Node metadata (244 unique proteins)
  - Columns: `node_id`, `gene_name`, `entrez_id`, `organism`

### 2. Enriched Data Files
- **`enriched_metadata.csv`** / **`enriched_metadata.json`** - Extended gene properties from NCBI
  - Additional columns: genomic location, gene type, functional summary, aliases, etc.

## Properties Available via Entrez ID

Using NCBI's Entrez Gene database, you can retrieve the following properties:

### Basic Information
- **Gene Symbol** - Official gene symbol (e.g., ACE2, ACTB)
- **Gene Name** - Full gene name
- **Aliases** - Alternative names and symbols
- **Gene Type** - protein-coding, ncRNA, pseudogene, etc.
- **Organism** - Species name

### Genomic Location
- **Chromosome** - Chromosomal location (e.g., X, 7, 17)
- **Map Location** - Cytogenetic band (e.g., Xp22.2, 7p22.1)
- **Genomic Coordinates** - Start and stop positions
- **RefSeq Accession** - Chromosome accession version

### Functional Annotations
- **Summary** - Detailed description of gene function
- **Gene Ontology (GO)** terms - Biological processes, molecular functions, cellular components
- **Pathways** - KEGG, Reactome pathway associations
- **Protein Domains** - Conserved domains and motifs

### Sequence Data
- **mRNA sequences** - RefSeq nucleotide accessions
- **Protein sequences** - RefSeq protein accessions
- **Exon/Intron structure** - Gene structure information

### Expression & Clinical
- **Tissue expression** - Expression patterns across tissues
- **Disease associations** - OMIM, ClinVar links
- **Phenotypes** - Associated phenotypes
- **Clinical significance** - Medical relevance

### External Database Links
- **PubMed** - Literature citations
- **UniProt** - Protein database
- **Ensembl** - Genome browser
- **OMIM** - Disease database
- **Model organism databases** - MGI, RGD, etc.

## Scripts

### `process_biogrid.py`
Processes raw BioGRID tab-delimited data into:
- Edge list CSV (network structure)
- Metadata CSV (node attributes)

**Usage:**
```bash
python3 process_biogrid.py
```

### `fetch_gene_properties.py`
Fetches detailed gene properties from NCBI using Entrez IDs.

**Key Features:**
- Demonstrates available Entrez fields
- Extracts common properties (location, type, summary)
- Respects NCBI rate limits (3 requests/second)
- Outputs CSV and JSON formats

**Usage:**
```bash
# Install dependencies
pip3 install -r requirements.txt

# Run demo (processes first 10 genes)
python3 fetch_gene_properties.py
```

**To process all genes**, edit the script and change:
```python
enriched_data = process_metadata_file(metadata_file, output_file, sample_size=None)
```

**Important:** Update your email in the script:
```python
Entrez.email = "your.email@example.com"  # Required by NCBI
```

## Data Source

- **BioGRID** - Biological General Repository for Interaction Datasets
- **NCBI Gene** - National Center for Biotechnology Information Gene database
- Original data: BIOGRID-GENE-4383915-5.0.251.tab3.txt

## Requirements

- Python 3.x
- Biopython (for NCBI Entrez queries)

## Example: ACE2 Gene Properties

```
Gene Symbol: ACE2
Gene Name: angiotensin converting enzyme 2
Organism: Homo sapiens
Chromosome: X
Map Location: Xp22.2
Genomic Position: chrX:15,518,196-15,607,210
Summary: The protein encoded by this gene belongs to the angiotensin-converting
         enzyme family... ACE2 is a functional receptor for the spike glycoprotein
         of SARS-CoV and SARS-CoV-2, the causative agent of COVID-19.
```

## API Resources

- **Biopython Entrez**: https://biopython.org/docs/dev/api/Bio.Entrez.html
- **NCBI E-utilities**: https://www.ncbi.nlm.nih.gov/books/NBK25501/
- **BioGRID**: https://thebiogrid.org/

## Notes

- The dataset appears to be SARS-CoV-2 related protein interactions
- 354 interactions between 244 unique proteins
- Mix of viral (SARS-CoV-2) and human proteins
- Rate limiting: NCBI allows 3 requests/second (10/second with API key)
