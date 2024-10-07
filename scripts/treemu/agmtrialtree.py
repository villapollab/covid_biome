import sys
from ete3 import Tree
import pandas as pd

def filter_tree_by_tax_ids(tree, tax_ids_set):
    """
    Recursively filter leaf nodes in the tree based on tax_id values.
    Remove leaf nodes whose tax_id is not in the specified tax_ids_set.
    """
    for leaf in tree.iter_leaves():
        if int(leaf.name) not in tax_ids_set:
            leaf.delete()

def rename_tree_tips_with_otu(tree, taxid_otu_map):
    """
    Recursively rename leaf nodes in the tree with corresponding OTU values.
    Replace each leaf node's tax_id name with the OTU value from the input map.
    """
    for leaf in tree.iter_leaves():
        current_tax_id = int(leaf.name)
        if current_tax_id in taxid_otu_map:
            new_tip_label = taxid_otu_map[current_tax_id]
            leaf.name = new_tip_label
        else:
            leaf.delete()

def main(tree_file, map_file, output_tree_file):
    # Read the reference tree
    reference_tree = Tree(tree_file)

    # Load the input map file into a DataFrame
    map_df = pd.read_csv(map_file)

    # Create a dictionary to map tax_id values to corresponding OTU values
    taxid_otu_map = dict(zip(map_df['tax_id'], map_df['emOTU']))

    # Set of tax_id values from the input map file
    tax_ids_in_map = set(map_df['tax_id'])

    # Filter tree to retain only leaf nodes with tax_id values present in the input map
    filter_tree_by_tax_ids(reference_tree, tax_ids_in_map)

    # Rename tree tip labels with corresponding OTU values from the input map
    rename_tree_tips_with_otu(reference_tree, taxid_otu_map)

    # Write the modified tree with filtered and updated tip labels to a newick file
    reference_tree.write(outfile=output_tree_file)

    print(f"Filtered and renamed tree saved to '{output_tree_file}' in Newick format.")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <input_tree_file> <input_map_file> <output_tree_file>")
        sys.exit(1)

    input_tree_file = sys.argv[1]
    input_map_file = sys.argv[2]
    output_tree_file = sys.argv[3]

    main(input_tree_file, input_map_file, output_tree_file)
