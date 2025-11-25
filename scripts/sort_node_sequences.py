import json

# Load the node to sequence mapping
with open('../processed_data/node_to_sequence.json', 'r') as f:
    node_to_sequence = json.load(f)

print(f"Loaded {len(node_to_sequence)} entries.")

# Sort by node_id as integer, but keep as strings
sorted_mapping = dict(sorted(node_to_sequence.items(), key=lambda item: int(item[0])))

# Save the sorted mapping
with open('../processed_data/node_to_sequence.json', 'w') as f:
    json.dump(sorted_mapping, f, indent=2)

print(f"Sorted and saved {len(sorted_mapping)} entries to node_to_sequence.json")
print(f"First 5 node IDs: {list(sorted_mapping.keys())[:5]}")
print(f"Last 5 node IDs: {list(sorted_mapping.keys())[-5:]}")
