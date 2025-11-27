#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_vcf_path> <reference_sample_name> <local_zenodo_ref_dir>"
    echo "Example: $0 /path/to/my_calls.vcf HG002 /home/user/my_zenodo_data/truvari_ref_vcf"
    exit 1
fi

# Input arguments
input_vcf_path=$1      # The path to the VCF file to be benchmarked
ref_name=$2            # The name of the reference sample (e.g., HG002)
zenodo_ref_base_dir=$3 # The local path to the downloaded /data_zenodo/truvari_ref_vcf/
default_refdist=100000

# Construct the sample-specific reference directory path
# For example, if $3 is /local/data/truvari_ref_vcf/, the script looks for HG002_* files there.
ref_dir="$zenodo_ref_base_dir"

# --- Setup Checks ---
# Check if the reference directory exists
if [ ! -d "$ref_dir" ]; then
    echo "No reference directory found at path: $ref_dir"
    echo "Please ensure the path to your local copy of '/data_zenodo/truvari_ref_vcf/' is correct."
    exit 1
fi

# --- Prepare Input VCF ---
echo "Preparing input VCF file..."

# Copy the input VCF file to the current directory for processing
cp "$input_vcf_path" .

# Extract the basename of the input VCF file
input_vcf_basename=$(basename "$input_vcf_path")

# Determine the final gzipped VCF name and process if needed
if [[ "$input_vcf_basename" != *.gz ]]; then
    # If the file is not gzipped, compress and remove the original
    echo "Compressing $input_vcf_basename with bgzip..."
    bgzip -c "$input_vcf_basename" > "${input_vcf_basename}.gz"
    rm "$input_vcf_basename" 
    input_vcf="${input_vcf_basename}.gz"
else
    # If it's already gzipped, use the existing file
    input_vcf="$input_vcf_basename"
fi

# Index the compressed VCF file with tabix
echo "Indexing $input_vcf with tabix..."
tabix -p vcf "$input_vcf"

# --- Run Truvari Benchmarking ---

echo "Searching for reference VCFs matching pattern: *.vcf.gz in $ref_dir"

# Find all reference VCFs that start with the reference sample name
# The find command ensures we only process files starting with the sample name
for ref_vcf in $(ls $ref_dir/*.vcf.gz); do
    
    base_name=$(basename "$ref_vcf" .vcf.gz)
    
    # Reset refdist to default for each loop iteration
    current_refdist=$default_refdist
    
    # Check for length stratification patterns in the base_name and set refdist accordingly
    if [[ "$base_name" == *50bp_10kb* ]]; then
        current_refdist=1000
        echo "INFO: Found '50bp_10kb' in $base_name. Setting refdist/chunksize to $current_refdist"
    elif [[ "$base_name" == *10kb_100kb* ]]; then
        current_refdist=10000
        echo "INFO: Found '10kb_100kb' in $base_name. Setting refdist/chunksize to $current_refdist"
    elif [[ "$base_name" == *100kb_1Mb* ]] || [[ "$base_name" == *1Mb_plus* ]]; then
        current_refdist=100000
        echo "INFO: Found '100kb_1Mb' or '1Mb_plus' in $base_name. Setting refdist/chunksize to $current_refdist"
    else
        echo "INFO: No specific length stratification found in $base_name. Using default refdist/chunksize of $current_refdist"
    fi

    # Define the output directory name based on the reference name and base name
    output_dir_name="${ref_name}_${base_name}"
    
    echo "Running truvari bench for $ref_vcf against $input_vcf..."

    # Run the first truvari bench with the dynamically set refdist and chunksize
    truvari bench --pctseq 0 --pick multi --chunksize "$current_refdist" --refdist "$current_refdist" --pctsize 0.3 -b "$ref_vcf" -c "$input_vcf" --sizemax 5400000 -o "./${output_dir_name}_pctseq0_sizemax_5.4mb"
    
    # Post-processing step
    output_path="./${output_dir_name}"
    if [ -d "$output_path" ]; then
        python cal_median_normalized_breakpoint_length_deviation.py "$output_path/tp-comp.vcf.gz"
        mv length_breakpoint_deviation.log "$output_path"
    else
        echo "WARNING: truvari bench failed to create output directory for $base_name. Skipping deviation calculation."
    fi

    python cal_tier_specific_stats.py ./ $ref_name --output ${ref_name}_tier_stats.json

    echo "Finished running truvari bench for $ref_vcf against $input_vcf."

done

echo "All truvari bench runs completed."
