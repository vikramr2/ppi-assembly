import json

def merge_dicts(dict1_path, dict2_path, output_path):
    """Merge two JSON dictionaries and save to output path."""
    with open(dict1_path, 'r') as f1:
        dict1 = json.load(f1)
    
    with open(dict2_path, 'r') as f2:
        dict2 = json.load(f2)
    
    # Merge dictionaries (dict2 entries will overwrite dict1 if keys overlap)
    merged_dict = {**dict1, **dict2}
    
    with open(output_path, 'w') as fout:
        json.dump(merged_dict, fout, indent=2)
    
    print(f"Merged dictionary saved to {output_path}")

dict1_path = '../processed_data/node_to_sequence.json'
dict2_path = '../processed_data/node_to_sequence_missing_uniprot.json'
output_path = '../processed_data/node_to_sequence.json'

merge_dicts(dict1_path, dict2_path, output_path)

# Remove the intermediate file if needed
import os
if os.path.exists(dict2_path):
    os.remove(dict2_path)
    print(f"Removed intermediate file: {dict2_path}")
