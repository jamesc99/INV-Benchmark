import re
import argparse

# Set up argument parsing
parser = argparse.ArgumentParser(description="Filter VCF by SVTYPE containing 'INV' and add SIMPLE_SV=INV.")
parser.add_argument("vcf_file", help="Path to the input VCF file.")
parser.add_argument("output_file", help="Path to the output VCF file.")

# Parse the arguments
args = parser.parse_args()

vcf_file = args.vcf_file
output_file = args.output_file

with open(vcf_file, 'r') as f_in, open(output_file, 'w') as f_out:
    header_added = False
    for line in f_in:
        if line.startswith("##") and not header_added:
            f_out.write(line)
            continue
        
        if line.startswith("#CHROM") and not header_added:
            # Add the new INFO header before the #CHROM line
            f_out.write('##INFO=<ID=SIMPLE_SV,Number=1,Type=String,Description="simple type of SV, only INV">\n')
            f_out.write(line)
            header_added = True
            continue

        if line.startswith("#"):
            f_out.write(line)
            continue
        
        columns = line.strip().split("\t")
        info_field = columns[7]  # INFO field is the 8th column
        if re.search(r'SVTYPE=.*INV', info_field):
            # Add SIMPLE_SV=INV to the INFO field if not already present
            if "SIMPLE_SV=INV" not in info_field:
                info_field += ";SIMPLE_SV=INV"
            columns[7] = info_field
            # Write the modified line to the new VCF file
            f_out.write("\t".join(columns) + "\n")
