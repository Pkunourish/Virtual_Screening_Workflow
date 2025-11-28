import argparse
import os
import re
import gzip
import pandas as pd
from typing import List, Optional, Tuple

def find_molecule_in_sdfgz(sdfgz_path: str, target_mol_name: str) -> Optional[List[str]]:
    """
    Search for a molecule by name in a gzipped SDF file.
    
    Args:
        sdfgz_path: Path to the .sdfgz file
        target_mol_name: Full molecule name to search for
        
    Returns:
        List of lines for the molecule block if found, None otherwise
    """
    try:
        with gzip.open(sdfgz_path, 'rt') as f_in:
            mol_block = []
            for line in f_in:
                mol_block.append(line)
                if line.strip() == '$$$$':
                    # Check if first line matches target molecule name
                    if mol_block[0].strip() == target_mol_name:
                        return mol_block
                    mol_block = []
    except Exception as e:
        raise RuntimeError(f"Error reading {sdfgz_path}: {str(e)}")
    return None

def construct_sdfgz_path(data_root: str, mol_name: str) -> str:
    """
    Construct the expected sdfgz file path from molecule name.
    
    Args:
        data_root: Root directory of data
        mol_name: Full molecule name
        
    Returns:
        Full path to the corresponding .sdfgz file
    
    Raises:
        ValueError: If molecule name format is invalid
    """
    pattern = r'^ligprep_(.*_split-\d+)_(\d+)\.sdf:\d+$'
    match = re.match(pattern, mol_name)
    if not match:
        raise ValueError(f"Invalid molecule name format: {mol_name}")
    
    base_part, split_num = match.groups()
    parent_dir = f"glide-dock_SP_{base_part}"
    sub_dir = f"glide-dock_SP_{base_part}_{split_num}"
    sdfgz_file = f"glide-dock_SP_{base_part}_{split_num}_lib.sdfgz"
    
    return os.path.join(data_root, parent_dir, sub_dir, sdfgz_file)

def process_csv(
    csv_path: str,
    data_root: str,
    output_sdf: str,
    mol_name_column: str = 'title'
) -> Tuple[int, List[str]]:
    """
    Process CSV file to extract molecules and save to SDF.
    
    Args:
        csv_path: Path to input CSV file
        data_root: Root directory of docking results
        output_sdf: Output SDF file path
        mol_name_column: Column name containing molecule identifiers
        
    Returns:
        Tuple containing (success_count, error_messages)
    """
    # Read CSV file
    df = pd.read_csv(csv_path)
    
    # Validate required column exists
    if mol_name_column not in df.columns:
        raise ValueError(f"Column '{mol_name_column}' not found in CSV")
    
    success_count = 0
    error_messages = []
    all_mol_blocks = []

    # Process each molecule entry
    for mol_name in df[mol_name_column]:
        try:
            # Construct SDFGZ path
            sdfgz_path = construct_sdfgz_path(data_root, mol_name)
            if not os.path.exists(sdfgz_path):
                raise FileNotFoundError(f"SDFGZ file not found: {sdfgz_path}")
            
            # Find molecule block
            mol_block = find_molecule_in_sdfgz(sdfgz_path, mol_name)
            if not mol_block:
                raise ValueError(f"Molecule not found in {sdfgz_path}")
            
            all_mol_blocks.append(mol_block)
            success_count += 1
            
        except Exception as e:
            error_messages.append(f"{mol_name}: {str(e)}")
            continue

    # Write collected molecules to SDF
    if all_mol_blocks:
        os.makedirs(os.path.dirname(output_sdf), exist_ok=True)
        with open(output_sdf, 'w') as f_out:
            for mol_block in all_mol_blocks:
                f_out.writelines(mol_block)

    return success_count, error_messages

def main():
    parser = argparse.ArgumentParser(
        description='Batch extract molecules from glide-dock results using CSV input'
    )
    parser.add_argument('--csv_path', required=True, help='Input CSV file path')
    parser.add_argument('--data_root', required=True, help='Root directory of docking results')
    parser.add_argument('--output_sdf', required=True, help='Output SDF file path')
    parser.add_argument('--mol_name_column', default='title',
                      help='CSV column containing molecule names (default: title)')
    args = parser.parse_args()

    try:
        success_count, errors = process_csv(
            args.csv_path,
            args.data_root,
            args.output_sdf,
            args.mol_name_column
        )
        
        print(f"Successfully extracted {success_count} molecules to {args.output_sdf}")
        
        if errors:
            print(f"\nEncountered {len(errors)} errors:")
            for error in errors[:5]:  # Show first 5 errors to avoid flooding
                print(f" - {error}")
            if len(errors) > 5:
                print(f" ...and {len(errors)-5} more errors")
                
    except Exception as e:
        print(f"Fatal error: {str(e)}")
        exit(1)

if __name__ == '__main__':
    main()