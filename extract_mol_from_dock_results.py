import argparse
import os
import re
import gzip
from typing import Optional

def find_molecule_in_sdfgz(sdfgz_path: str, target_mol_name: str) -> Optional[list]:
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
    # Extract components from molecule name using regex
    pattern = r'^ligprep_(.*_split-\d+)_(\d+)\.sdf:\d+$'
    match = re.match(pattern, mol_name)
    if not match:
        raise ValueError(f"Invalid molecule name format: {mol_name}")
    
    base_part, split_num = match.groups()
    parent_dir = f"glide-dock_SP_{base_part}"
    sub_dir = f"glide-dock_SP_{base_part}_{split_num}"
    sdfgz_file = f"glide-dock_SP_{base_part}_{split_num}_lib.sdfgz"
    
    return os.path.join(data_root, parent_dir, sub_dir, sdfgz_file)

def main():
    parser = argparse.ArgumentParser(description='Extract specific molecule from glide-dock results')
    parser.add_argument('--data_root', required=True, help='Root directory containing glide-dock results')
    parser.add_argument('--mol_name', required=True, help='Full molecule name to extract (e.g., "ligprep_chemdiv-stock1_04_414266_split-10_5.sdf:5784")')
    parser.add_argument('--output_dir', required=True, help='Directory to save extracted molecule')
    args = parser.parse_args()

    # Create output directory if needed
    os.makedirs(args.output_dir, exist_ok=True)

    try:
        # Construct expected sdfgz file path
        sdfgz_path = construct_sdfgz_path(args.data_root, args.mol_name)
        if not os.path.exists(sdfgz_path):
            raise FileNotFoundError(f"SDFGZ file not found: {sdfgz_path}")

        # Search for molecule
        mol_block = find_molecule_in_sdfgz(sdfgz_path, args.mol_name)
        if not mol_block:
            raise ValueError(f"Molecule {args.mol_name} not found in {sdfgz_path}")

        # Create safe filename and save
        safe_name = args.mol_name.replace(':', '_').replace('/', '_')
        output_path = os.path.join(args.output_dir, f"{safe_name}.sdf")
        with open(output_path, 'w') as f_out:
            f_out.writelines(mol_block)
        
        print(f"Successfully saved molecule to: {output_path}")

    except Exception as e:
        print(f"Error: {str(e)}")
        exit(1)

if __name__ == '__main__':
    main()