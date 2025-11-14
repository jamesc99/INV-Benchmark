#!/bin/bash
#SBATCH --job-name=truvari_bench_all
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --time=3:00:00
#SBATCH -A proj-fs0002

# Input arguments
input_vcf=$1
ref_name=$2
default_refdist=100000

# Reference directory based on the input reference name
ref_dir=/users/u250191/ryan_scratch_ln/benchmark_inv/reference_Cell_paper/assembly_based_category/truvari_ref/round1_try2simple/$ref_name

# Check if the corresponding reference directory exists
if [ ! -d "$ref_dir" ]; then
    echo "No reference directory found for reference: $ref_name. Exiting."
    exit 1
fi

# Source the bashrc to get the conda environment
. /users/u250191/.bashrc
mamba activate truvari

# Copy the input VCF file to the current directory
cp $input_vcf .

# Extract the basename of the input VCF file
input_vcf_basename=$(basename $input_vcf)

# Compress the VCF file with bgzip
bgzip -c $input_vcf_basename > ${input_vcf_basename}.gz

# Index the compressed VCF file with tabix
tabix -p vcf ${input_vcf_basename}.gz

# Update the input VCF file to the compressed version for further processing
input_vcf="${input_vcf_basename}.gz"

# Loop through each VCF file in the reference folder
for ref_vcf in $(ls $ref_dir/*.vcf.gz); do
    base_name=$(basename $ref_vcf .vcf.gz)

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
    truvari bench --pctseq 0 --pick multi --chunksize $current_refdist --refdist $current_refdist --pctsize 0.3 -b $ref_vcf -c $input_vcf --sizemax 5400000 -o ./${output_dir_name}_pctseq0_sizemax_5.4mb
    python /users/u250191/ryan_scratch_ln/scripts/benchmark_inv/truvari/cal_median_normalized_breakpoint_length_deviation.py ./${output_dir_name}_pctseq0_sizemax_5.4mb/tp-comp.vcf.gz
    mv length_breakpoint_deviation.log ./${output_dir_name}_pctseq0_sizemax_5.4mb

    # Run the second truvari bench with --pctseq 0
    #truvari bench --pctseq 0 --pick multi --refdist 1000 --pctsize 0.3 -b $ref_vcf -c $input_vcf -o ./${output_dir_name}_pctseq0_sizemax_50kb
    #python /users/u250191/ryan_scratch_ln/scripts/benchmark_inv/truvari/cal_breakpoint_length_deviation.py ./${output_dir_name}_pctseq0_sizemax_50kb/tp-comp.vcf.gz
    #mv length_breakpoint_deviation.log ./${output_dir_name}_pctseq0_sizemax_50kb

    echo "Finished running truvari bench for $ref_vcf against $input_vcf."
done

echo "All truvari bench runs completed."
