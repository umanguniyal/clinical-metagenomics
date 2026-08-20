import sys

def convert_kreport(input_file, output_file):
    lineage = []
    with open(input_file) as f, open(output_file, 'w') as out:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 6:
                continue
            
            # parts[2] is the number of reads assigned directly to this taxon
            count = int(parts[2])
            name = parts[5]
            
            # calculate depth by counting leading spaces (Kraken2 uses 2 spaces per level)
            leading_spaces = len(name) - len(name.lstrip(' '))
            depth = leading_spaces // 2
            
            name = name.strip()
            
            # Adjust the lineage based on current depth
            if depth < len(lineage):
                lineage = lineage[:depth]
            lineage.append(name)
            
            # Only write to krona if reads are assigned directly here
            if count > 0:
                out.write(f"{count}\t" + "\t".join(lineage) + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python kreport2text.py <input_kreport> <output_text>")
        sys.exit(1)
    convert_kreport(sys.argv[1], sys.argv[2])
