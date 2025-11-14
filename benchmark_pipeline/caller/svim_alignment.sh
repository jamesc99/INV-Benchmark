#!/bin/bash
#SBATCH --job-name=svim
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=72:00:00
#SBATCH --partition=medium

# --- Environment Setup ---
mamba activate svim_env
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
ref38="human-grch38.fasta"
WORK_DIR=$(pwd)

if [ -z "$1" ]; then
    echo "Error: No input BAM file provided."
    echo "Usage: sbatch $0 /path/to/your/file.bam"
    exit 1
fi
input_bam=$1

# --- Automatic Naming and Directory Setup ---
if [ -z "$SLURM_JOB_ID" ]; then
    echo "Error: This script must be run as a Slurm job." >&2
    exit 1
fi
job_id=$SLURM_JOB_ID

# Create a unique output directory for this job to prevent conflicts
OUT_DIR="${WORK_DIR}/svim_out_${job_id}"
mkdir -p "$OUT_DIR"
echo "Using Slurm Job ID ${job_id}. Output will be in: ${OUT_DIR}"


# --- Main Execution & Time Logging ---
TIME_LOG_TMP=$(mktemp)
echo "Starting SVIM analysis..."
# The time command's stderr is redirected to the temporary file.
# SVIM will write its output to the unique OUT_DIR.
{ time svim alignment ${OUT_DIR} ${input_bam} ${ref38}; } 2> "$TIME_LOG_TMP"
echo "SVIM analysis complete."


# --- Process and Record Running Time ---
echo "Processing and recording running time..."
real_time=$(grep '^real' "$TIME_LOG_TMP" | awk '{print $2}')
user_time_str=$(grep '^user' "$TIME_LOG_TMP" | awk '{print $2}')
sys_time_str=$(grep '^sys' "$TIME_LOG_TMP" | awk '{print $2}')

# Convert user and sys times to seconds
user_seconds=$(time_to_seconds "$user_time_str")
sys_seconds=$(time_to_seconds "$sys_time_str")

# Calculate the total CPU time
total_cpu_seconds=$(echo "$user_seconds + $sys_seconds" | bc)

# Append the formatted time information to the main running_time.log file.
(
  echo "job_id: ${job_id}"
  echo "total_cpu_time: ${total_cpu_seconds}s (${user_seconds}s user + ${sys_seconds}s sys)"
  echo "wall_clock_time: ${real_time}"
  echo "---"
) >> running_time.log

# Clean up the temporary file.
rm "$TIME_LOG_TMP"
echo "Running time successfully logged to running_time.log"


# --- Post-processing ---
# All post-processing now uses files from the unique output directory
echo "Filtering VCF with bcftools..."
bcftools view -i 'QUAL >= 10' ${OUT_DIR}/variants.vcf > ${OUT_DIR}/variants_qual_greater10.vcf

echo "Filtering for inversions..."
grep -E '^#|SVTYPE=INV' ${OUT_DIR}/variants_qual_greater10.vcf > ${OUT_DIR}/filtered_qual_10_inv.vcf
grep -v '#' ${OUT_DIR}/filtered_qual_10_inv.vcf > ${OUT_DIR}/filtered_qual_10_inv_woheader.vcf
echo "All steps complete."

