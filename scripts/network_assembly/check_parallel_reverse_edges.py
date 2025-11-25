import csv

# Read the edge list
edges = []
with open('brca_ppi_edgelist.csv', 'r') as f:
    reader = csv.reader(f)
    header = next(reader)  # Skip header
    for row in reader:
        source, target = int(row[0]), int(row[1])
        edges.append((source, target))

print(f"Total edges: {len(edges)}")

# Remove duplicates and reverse edges
# Normalize edges so that source <= target (for undirected graph)
normalized_edges = set()
for source, target in edges:
    # Store as (min, max) to treat as undirected
    normalized = (min(source, target), max(source, target))
    normalized_edges.add(normalized)

print(f"Unique undirected edges after removing duplicates and reverse edges: {len(normalized_edges)}")
print(f"Removed: {len(edges) - len(normalized_edges)} edges")

# Write cleaned edge list
output_file = 'brca_ppi_edgelist_cleaned.csv'
with open(output_file, 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['source', 'target'])
    for source, target in sorted(normalized_edges):
        writer.writerow([source, target])

print(f"\nCleaned edge list written to: {output_file}")
print(f"Format: Each edge appears only once with source <= target")
