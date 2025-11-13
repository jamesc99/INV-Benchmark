import sys
import re
import os

def append_end_to_vcf(vcf_file):
    # Extract the base name of the input file (without extension)
    base_name = os.path.basename(vcf_file).replace(".vcf", "")
    
    # File names for the output files
    output_with_end = f"add_end_{base_name}.vcf"
    output_without_header = f"add_end_{base_name}_woheader.vcf"
    
    with open(vcf_file, 'r') as infile, open(output_with_end, 'w') as outfile_with_end, open(output_without_header, 'w') as outfile_woheader:
        for line in infile:
            # Keep header lines unchanged in the first file, but ignore them in the second file
            if line.startswith("#"):
                outfile_with_end.write(line)
                continue
            
            # Split the line by tabs
            fields = re.split(r'\s+', line.strip())
            if len(fields) < 8:
                outfile_with_end.write(line)
                continue

            chr_val = fields[0]
            pos = int(fields[1])
            info = fields[7]

            # Extract SVLEN from the INFO field
            svlen_match = re.search(r'SVLEN=(\d+)', info)
            svlen = int(svlen_match.group(1)) if svlen_match else 0
            end = pos + svlen

            # Append END field to INFO if not already present
            if 'END=' not in info:
                info += f";END={end}"

            # Update the INFO field and write the modified line to the first output file
            fields[7] = info
            outfile_with_end.write("\t".join(fields) + "\n")

            # Write the same line to the second output file without headers
            outfile_woheader.write("\t".join(fields) + "\n")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python script.py <input_vcf>")
        sys.exit(1)

    vcf_file = sys.argv[1]
    append_end_to_vcf(vcf_file)