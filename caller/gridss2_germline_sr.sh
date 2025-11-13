#!/bin/bash
#SBATCH --job-name=gridss2_germline
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=48:00:00
#SBATCH --partition=medium

# --- Environment Setup ---
mamba activate gridss-2.13.2
set -eox # Exit on error, print commands

# --- Helper Function for Time Conversion ---
# Converts time format (e.g., 1m2.345s) into total seconds (e.g., 62.345)
time_to_seconds() {
    local time_str=$1
    local minutes=0
    local seconds=0

    if [[ "$time_str" == *m* ]]; then
        minutes=$(echo "$time_str" | cut -d'm' -f1)
        seconds_part=$(echo "$time_str" | cut -d'm' -f2)
        seconds=$(echo "$seconds_part" | sed 's/s//')
    else
        seconds=$(echo "$time_str" | sed 's/s//')
    fi
    # Use bc for floating-point arithmetic
    echo "$minutes * 60 + $seconds" | bc
}

# --- Configuration & Argument Handling ---
exclude_bed_hg38="hg38_exclude_list_ENCFF356LFX.bed"
reference="human-grch38.fasta"
working_dir=$(pwd)

if [ -z "$1" ]; then
    echo "Error: No input BAM file provided."
    echo "Usage: sbatch $0 /path/to/your/file.bam"
    exit 1
fi
input_bam=$1

# --- Automatic Naming ---
if [ -z "$SLURM_JOB_ID" ]; then
    echo "Error: This script must be run as a Slurm job." >&2
    exit 1
fi
job_id=$SLURM_JOB_ID
output_vcf="${job_id}.gridss.vcf"
echo "Using Slurm Job ID ${job_id}. Output VCF will be: ${output_vcf}"

# --- Main Execution & Time Logging ---
TIME_LOG_TMP=$(mktemp)
echo "Starting GRIDSS2 analysis..."
# The time command's stderr is redirected to the temporary file.
{ time gridss \
    -r ${reference} \
    -t 8 \
    -b ${exclude_bed_hg38} \
    -o ${output_vcf} \
    ${input_bam}; } 2> "$TIME_LOG_TMP"
echo "GRIDSS2 analysis complete."


# --- Process and Record Running Time ---
echo "Processing and recording running time..."
real_time=$(grep '^real' "$TIME_LOG_TMP" | awk '{print $2}')
user_time_str=$(grep '^user' "$TIME_LOG_TMP" | awk '{print $2}')
sys_time_str=$(grep '^sys' "$TIME_LOG_TMP" | awk '{print $2}')

user_seconds=$(time_to_seconds "$user_time_str")
sys_seconds=$(time_to_seconds "$sys_time_str")
total_cpu_seconds=$(echo "$user_seconds + $sys_seconds" | bc)

(
  echo "job_id: ${job_id}"
  echo "total_cpu_time: ${total_cpu_seconds}s (${user_seconds}s user + ${sys_seconds}s sys)"
  echo "wall_clock_time: ${real_time}"
  echo "---"
) >> running_time.log

rm "$TIME_LOG_TMP"
echo "Running time successfully logged to running_time.log"


# --- Post-processing ---
mamba deactivate
mamba activate R_env_4.4.0

echo "Running downstream R script for inversion analysis..."
# Run the external R script on the uniquely named output VCF
Rscript /users/u250191/ryan_scratch_ln/benchmark_inv/caller/gridss2/downstream_inv.R "${output_vcf}" "hg38" "${working_dir}"
echo "All steps complete."

