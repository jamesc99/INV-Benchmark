#!/usr/bin/env python3
import sys
import os

def vcf_to_bed(vcf_file):
    # Open the VCF file for reading
    with open(vcf_file, 'r') as vcf:
        # Open the BED file for writing (use the same name with .bed extension)
        bed_file = vcf_file.replace('.vcf', '.bed')
        with open(bed_file, 'w') as bed:
            for line in vcf:
                # Skip header lines in the VCF file
                if line.startswith("#"):
                    continue
                # Split the VCF line into fields
                fields = line.strip().split('\t')
                chrom = fields[0]  # CHROM column
                start = fields[1]  # POS column
                info = fields[7]   # INFO column
                
                # Extract END value from the INFO field
                info_fields = info.split(';')
                end = None
                for field in info_fields:
                    if field.startswith("END="):
                        end = field.split('=')[1]
                        break
                
                if end:
                    # Write the BED line in format: chrom, start-1 (0-based), end
                    bed.write(f"{chrom}\t{int(start) - 1}\t{end}\n")
                else:
                    print(f"Warning: No END found for line: {line.strip()}")

def process_subdirectories():
    # Walk through current directory and its subdirectories
    for root, dirs, files in os.walk('.'):
        # Target only inv_severus.vcf in each subdirectory
        if 'inv_severus.vcf' in files:
            vcf_file_path = os.path.join(root, 'add_end_inv_severus.vcf')
            print(f"Processing file: {vcf_file_path}")
            vcf_to_bed(vcf_file_path)
            print(f"Converted {vcf_file_path} to BED format")

if __name__ == "__main__":
    process_subdirectories()
