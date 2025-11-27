#!/bin/bash

# Check if the correct number of arguments is provided
# We need paths to the python scripts to be safe, so I added arguments 4 and 5
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_vcf_path> <reference_sample_name> <local_zenodo_ref_dir>"
    echo "Example: $0 calls.vcf HG002 /data/ref_vcf"
    exit 1
fi

# Input arguments
input_vcf_path=$1       
ref_name=$2             
zenodo_ref_base_dir=$3  
default_refdist=100000

ref_dir="$zenodo_ref_base_dir"

# --- Setup Checks ---
if [ ! -d "$ref_dir" ]; then
    echo "Error: Reference directory not found at: $ref_dir"
    exit 1
fi

# --- Prepare Input VCF ---
echo "Preparing input VCF file..."
cp "$input_vcf_path" .
input_vcf_basename=$(basename "$input_vcf_path")

if [[ "$input_vcf_basename" != *.gz ]]; then
    echo "Compressing $input_vcf_basename with bgzip..."
    bgzip -c "$input_vcf_basename" > "${input_vcf_basename}.gz"
    rm "$input_vcf_basename" 
    input_vcf="${input_vcf_basename}.gz"
else
    input_vcf="$input_vcf_basename"
fi

echo "Indexing $input_vcf with tabix..."
tabix -p vcf "$input_vcf"

# --- Run Truvari Benchmarking ---

echo "Searching for reference VCFs matching pattern: ${ref_name}*.vcf.gz in $ref_dir"

for ref_vcf in $(ls $ref_dir/*.vcf.gz); do
    base_name=$(basename "$ref_vcf" .vcf.gz)
    current_refdist=$default_refdist
    if [[ "$base_name" == *50bp_10kb* ]]; then
        current_refdist=1000
        echo "INFO: Found '50bp_10kb' in $base_name. Setting refdist to $current_refdist"
    elif [[ "$base_name" == *10kb_100kb* ]]; then
        current_refdist=10000
        echo "INFO: Found '10kb_100kb' in $base_name. Setting refdist to $current_refdist"
    elif [[ "$base_name" == *100kb_1Mb* ]] || [[ "$base_name" == *1Mb_plus* ]]; then
        current_refdist=100000
        echo "INFO: Found '100kb_1Mb' in $base_name. Setting refdist to $current_refdist"
    else
        echo "INFO: Using default refdist $current_refdist"
    fi
    output_folder_name="${ref_name}_${base_name}"
    
    echo "Running truvari bench for $ref_vcf..."

    truvari bench --pctseq 0 --pick multi \
        --chunksize "$current_refdist" \
        --refdist "$current_refdist" \
        --pctsize 0.3 \
        -b "$ref_vcf" \
        -c "$input_vcf" \
        --sizemax 5400000 \
        -o "./$output_folder_name"
    
    # FIX 3: Check the CORRECT directory name
    output_path="./$output_folder_name"
    
    if [ -d "$output_path" ]; then
        # Run deviation calculation
        python benchmark_pipeline/truvari/cal_median_normalized_breakpoint_length_deviation.py "$output_path/tp-comp.vcf.gz"
        mv length_breakpoint_deviation.log "$output_path"
    else
        echo "WARNING: Truvari output directory not found at $output_path"
    fi

    echo "Finished running truvari bench for $base_name."

done

echo "All Truvari runs completed. Calculating Tier Specific Stats..."

python benchmark_pipeline/truvari/cal_tier_specific_stats.py ./ "$ref_name" --output "${ref_name}_tier_stats.json"

echo "Pipeline Finished."
