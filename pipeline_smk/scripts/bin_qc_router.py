import argparse
import os
import shutil

def get_kingdom(bin_tax_file):
    # A placeholder for kingdom detection.
    # Reads bin taxonomy from GTDB-Tk or Kraken2 output,
    # returns 'bacteria', 'archaea', 'viral', or 'eukarya'
    # Fallback to bacteria if unknown
    
    if not os.path.exists(bin_tax_file):
        return 'bacteria'
        
    with open(bin_tax_file, 'r') as f:
        content = f.read().lower()
        if 'virus' in content or 'viral' in content or 'viricota' in content:
            return 'viral'
        if 'eukary' in content or 'fungi' in content:
            return 'eukarya'
        if 'archaea' in content:
            return 'archaea'
            
    return 'bacteria'

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin_dir", required=True)
    parser.add_argument("--bin_tax_dir", required=True)
    parser.add_argument("--out_checkm2", required=True)
    parser.add_argument("--out_checkv", required=True)
    parser.add_argument("--out_busco", required=True)
    
    args = parser.parse_args()
    
    os.makedirs(args.out_checkm2, exist_ok=True)
    os.makedirs(args.out_checkv, exist_ok=True)
    os.makedirs(args.out_busco, exist_ok=True)
    
    for bin_file in os.listdir(args.bin_dir):
        if not bin_file.endswith('.fa'):
            continue
            
        bin_path = os.path.join(args.bin_dir, bin_file)
        bin_name = bin_file.replace('.fa', '')
        
        # Determine kingdom
        tax_file = os.path.join(args.bin_tax_dir, f"{bin_name}.best_hit.txt")
        kingdom = get_kingdom(tax_file)
        
        # Route
        if kingdom in ['bacteria', 'archaea']:
            shutil.copy(bin_path, os.path.join(args.out_checkm2, bin_file))
        elif kingdom == 'viral':
            shutil.copy(bin_path, os.path.join(args.out_checkv, bin_file))
        elif kingdom == 'eukarya':
            shutil.copy(bin_path, os.path.join(args.out_busco, bin_file))

if __name__ == "__main__":
    main()
