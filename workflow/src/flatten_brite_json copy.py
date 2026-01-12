#!/usr/bin/env python3
"""
Transform nested KEGG BRITE JSON to a flat table.
Each leaf term becomes a row with its full hierarchical path.
"""
# %%
import json
import re
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
import sys

# %%
def extract_leaf_paths(node, path=None):
    """
    Recursively traverse the JSON tree and extract all leaf nodes with their paths.
    
    Parameters:
    -----------
    node : dict
        Current node in the tree
    path : list
        Current path from root to this node
        
    Yields:
    -------
    list : Full path from root to leaf node
    """
    if path is None:
        path = []
    
    # Get the current node's name
    current_name = node.get("name", "")
    current_path = path + [current_name]
    
    # Check if this is a leaf node (no children or empty children)
    children = node.get("children", [])
    
    if not children:
        # This is a leaf node
        yield current_path
    else:
        # Recursively process children
        for child in children:
            yield from extract_leaf_paths(child, current_path)


def flatten_brite_json(json_file, output_file=None, max_depth=None):
    """
    Convert nested BRITE JSON to a flat table.
    
    Parameters:
    -----------
    json_file : str or Path
        Path to input JSON file
    output_file : str or Path, optional
        Path to output TSV file. If None, uses same name with .tsv extension
    max_depth : int, optional
        Maximum depth to use for columns. If None, uses the maximum depth found.
    
    Returns:
    --------
    pd.DataFrame : Flattened table
    """
    # Load JSON
    with open(json_file, 'r') as f:
        data = json.load(f)
    
    # Extract all leaf paths
    leaf_paths = list(extract_leaf_paths(data))
    
    # Determine maximum depth
    if max_depth is None:
        max_depth = max(len(path) for path in leaf_paths)
    
    # Create column names
    columns = ["root"] + [f"level_{i}" for i in range(1, max_depth)]
    
    # Build rows
    rows = []
    for path in leaf_paths:
        # Pad path to max_depth if needed
        name_pattern = re.compile(r'^\d+ ([^;]+);')
        elements_before_name = []
        name_elements = ""
        for element in path:
            match = name_pattern.match(element)
            if match:
                name_elements = element
            else:
                elements_before_name.append(element)        
        padded_path = elements_before_name + [""] * (max_depth - len(path)) + [name_elements]
        rows.append(padded_path[:max_depth])
    
    # Create DataFrame
    df = pd.DataFrame(rows, columns=columns)
    
    # Split the last level (leaf node) by tab character if present
    # Find the rightmost non-empty column (the leaf level)
    leaf_col = None
    for col in reversed(df.columns):
        if df[col].str.len().sum() > 0:
            leaf_col = col
            break
    
    if leaf_col:
        # Split by tab and expand into separate columns
        split_df = df[leaf_col].str.split('\t', expand=True)
        
        # Determine how many split columns we have
        num_splits = split_df.shape[1]
        
        if num_splits > 1:
            # Rename the split columns
            split_cols = [f"{leaf_col}_part_{i+1}" for i in range(num_splits)]
            split_df.columns = split_cols
            
            # Extract gene name from the first part (format: "code gene_name; description")
            first_col = split_cols[0]
            if first_col in split_df.columns:
                # Extract gene name: text between first space and semicolon
                gene_names = split_df[first_col].str.extract(r'^\d+\s+([^;]+);', expand=False)
                # Remove "SPOM_" prefix if present
                gene_names = gene_names.str.replace(r'^SPOM_', '', regex=True)
                split_df.insert(1, 'gene_name', gene_names)
            
            # Remove the original leaf column and add the split columns
            df = df.drop(columns=[leaf_col])
            df = pd.concat([df, split_df], axis=1)
            
            print(f"Split leaf terms into {num_splits} columns with gene names extracted")
    
    # Reorder columns to put gene_name first
    if 'gene_name' in df.columns:
        cols = ['gene_name'] + [col for col in df.columns if col != 'gene_name']
        df = df[cols]
    
    col_dict = {i: col for i, col in enumerate(df.columns.tolist())}
    max_level = max(col_dict.keys())
    col_dict[max_level] = "level_8"
    col_dict[max_level - 1] = "level_7"
    for i in range(1, max_level - 1):
        col_dict[i] = f"level_{i}"
    df.columns = list(col_dict.values())
    for i in range(1,9):
        if f"level_{i}" not in df.columns:
            df[f"level_{i}"] = np.nan
    


    df = df.sort_index(axis=1)

    # Save to file if specified
    if output_file is None:
        output_file = Path(json_file).with_suffix(".tsv")
    
    df.to_csv(output_file, sep="\t", index=False, header=False)
    print(f"Wrote {len(df)} leaf terms to {output_file}")
    print(f"Hierarchy depth: {max_depth} levels")
    
    return df


def process_directory(input_dir, output_dir=None, max_depth=None):
    """
    Process all JSON files in a directory.
    
    Parameters:
    -----------
    input_dir : Path
        Input directory containing JSON files
    output_dir : Path, optional
        Output directory for TSV files. If None, uses input directory.
    max_depth : int, optional
        Maximum hierarchy depth for columns
        
    Returns:
    --------
    list : List of output file paths
    """
    input_path = Path(input_dir)
    
    if not input_path.is_dir():
        raise ValueError(f"Input path is not a directory: {input_dir}")
    
    # Find all JSON files
    json_files = list(input_path.glob("*.json"))
    
    if not json_files:
        print(f"No JSON files found in {input_dir}")
        return []
    
    print(f"Found {len(json_files)} JSON file(s) to process")
    print()
    
    # Set up output directory
    if output_dir is None:
        output_path = input_path
    else:
        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)
    
    output_files = []
    
    for i, json_file in enumerate(json_files, 1):
        print(f"[{i}/{len(json_files)}] Processing {json_file.name}...")
        
        output_file = output_path / json_file.with_suffix(".tsv").name
        
        try:
            df = flatten_brite_json(
                json_file,
                output_file=output_file,
                max_depth=max_depth
            )
            output_files.append(output_file)
            print()
        except Exception as e:
            print(f"Error processing {json_file.name}: {e}")
            print()
            continue
    
    print(f"\nSuccessfully processed {len(output_files)}/{len(json_files)} files")
    return output_files


def main():
    parser = argparse.ArgumentParser(
        description="Flatten nested KEGG BRITE JSON to a table",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process single file
  python flatten_brite_json.py brite_json/spo00001.json
  
  # Process single file with custom output
  python flatten_brite_json.py brite_json/spo00001.json -o output.tsv
  
  # Process all JSON files in a directory
  python flatten_brite_json.py brite_json/
  
  # Process directory with custom output directory
  python flatten_brite_json.py brite_json/ -o output_dir/
  
  # Limit hierarchy depth
  python flatten_brite_json.py brite_json/spo00001.json --max-depth 5
        """
    )
    
    parser.add_argument(
        "input",
        type=str,
        help="Input BRITE JSON file or directory containing JSON files"
    )
    
    parser.add_argument(
        "-o", "--output",
        type=str,
        default=None,
        help="Output TSV file or directory (default: same location as input)"
    )
    
    parser.add_argument(
        "--max-depth",
        type=int,
        default=None,
        help="Maximum hierarchy depth for columns (default: auto-detect)"
    )
    
    args = parser.parse_args()
    
    input_path = Path(args.input)
    
    # Check if input exists
    if not input_path.exists():
        print(f"Error: Input path does not exist: {args.input}", file=sys.stderr)
        sys.exit(1)
    
    # Process based on input type
    if input_path.is_dir():
        # Process directory
        output_files = process_directory(
            input_path,
            output_dir=args.output,
            max_depth=args.max_depth
        )
    else:
        # Process single file
        if args.output and Path(args.output).is_dir():
            # If output is a directory, construct output filename
            output_file = Path(args.output) / input_path.with_suffix(".tsv").name
            Path(args.output).mkdir(parents=True, exist_ok=True)
        else:
            output_file = args.output
        
        df = flatten_brite_json(
            input_path,
            output_file=output_file,
            max_depth=args.max_depth
        )
        
        # Display sample
        print("\nFirst 5 rows:")
        print(df.head().to_string())
        
        print(f"\nLast 5 rows:")
        print(df.tail().to_string())

# %%
if __name__ == "__main__":
    main()
