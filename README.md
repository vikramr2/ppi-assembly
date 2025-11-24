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
