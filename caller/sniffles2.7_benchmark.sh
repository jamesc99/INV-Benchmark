#!/bin/bash
#SBATCH --job-name=sniffle2
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=72:00:00
#SBATCH --partition=medium

mamba activate sniffles2_270
set -eox

# --- Helper Function for Time Conversion ---
# Converts time format (e.g., 1m2.345s) into total seconds (e.g., 62.345)
time_to_seconds() {
    local time_str=$1
    local minutes=0
    local seconds=0

    # Check if 'm' (minutes) is in the string
    if [[ "$time_str" == *m* ]]; then
        minutes=$(echo "$time_str" | cut -d'm' -f1)
        seconds_part=$(echo "$time_str" | cut -d'm' -f2)
        seconds=$(echo "$seconds_part" | sed 's/s//')
    else
        # No minutes, just seconds
        seconds=$(echo "$time_str" | sed 's/s//')
    fi

    # Use bc for floating-point arithmetic to calculate total seconds
    echo "$minutes * 60 + $seconds" | bc
}

# --- Reference Genome Paths ---
ref38="human-grch38.fasta"

# --- Argument and Name Handling ---
if [ "$#" -ne 1 ]; then
    echo "Illegal number of parameters. This script requires exactly one argument."
    echo "Usage: sbatch $0 /path/to/your/file.bam"
    exit 111
fi

bam_file=$1
threads=8

if [ -z "$SLURM_JOB_ID" ]; then
    echo "Error: This script must be run as a Slurm job." >&2
    echo "The SLURM_JOB_ID environment variable is not set." >&2
    exit 1
fi
name=$SLURM_JOB_ID
echo "Using Slurm Job ID as sample name: $name"

# --- Main Execution & Time Logging ---
TIME_LOG_TMP=$(mktemp)
echo "Starting Sniffles2 analysis..."
{ time sniffles --input $bam_file --vcf $name.vcf --sample-id ${name} --reference $ref38 --threads $threads > ${bam_file/.bam/.log} 2>&1; } 2> "$TIME_LOG_TMP"
echo "Sniffles2 analysis complete. Main output: $name.vcf"


# --- Process and Record Running Time ---
echo "Processing and recording running time..."
# Extract the 'real', 'user', and 'sys' time strings from the temp file.
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
  echo "job_id: ${name}"
  echo "total_cpu_time: ${total_cpu_seconds}s (${user_seconds}s user + ${sys_seconds}s sys)"
  echo "wall_clock_time: ${real_time}"
  echo "---"
) >> running_time.log

# Clean up the temporary file.
rm "$TIME_LOG_TMP"
echo "Running time successfully logged to running_time.log"


# --- Post-processing ---
echo "Filtering for inversions..."
grep -E '^#|SVTYPE=INV' $name.vcf > inv_sniffle2_${name}.vcf
grep -v '#' inv_sniffle2_${name}.vcf > inv_sniffle2_${name}_woheader.vcf
echo "Inversion filtering complete. Final files: inv_sniffle2_${name}.vcf, inv_sniffle2_${name}_woheader.vcf"

