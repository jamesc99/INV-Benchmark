#!/bin/bash
#SBATCH --job-name=cutesv2
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=72:00:00
#SBATCH --partition=medium

# --- Environment Setup ---
mamba activate cutesv2
set -eox # Exit on error, print commands

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


# --- Configuration ---
ref38="human-grch38.fasta"
WORK_DIR=$(pwd)

# --- Argument and Name Handling ---
# Check that exactly two arguments are provided: BAM file and data type.
if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters. This script requires exactly two arguments."
    echo "Usage: sbatch $0 /path/to/your/file.bam <pacbio|ont>"
    exit 1
fi

input_bam=$1
datatype=$2

# Automatically assign the Slurm Job ID to 'outputname'.
if [ -z "$SLURM_JOB_ID" ]; then
    echo "Error: This script must be run as a Slurm job." >&2
    echo "The SLURM_JOB_ID environment variable is not set." >&2
    exit 1
fi
outputname=$SLURM_JOB_ID
echo "Using Slurm Job ID as output name: $outputname"


# --- Set cuteSV Parameters ---
if [ "$datatype" == "pacbio" ]; then
    MAX_CLUSTER_BIAS_INS=1000
    DIFF_RATIO_MERGING_INS=0.9
    MAX_CLUSTER_BIAS_DEL=1000
    DIFF_RATIO_MERGING_DEL=0.5
elif [ "$datatype" == "ont" ]; then
    MAX_CLUSTER_BIAS_INS=100
    DIFF_RATIO_MERGING_INS=0.3
    MAX_CLUSTER_BIAS_DEL=100
    DIFF_RATIO_MERGING_DEL=0.3
else
    echo "Invalid data type specified. Please use 'pacbio' or 'ont'."
    exit 1
fi


# --- Main Execution & Time Logging ---
TIME_LOG_TMP=$(mktemp)
echo "Starting cuteSV analysis..."
{ time cuteSV \
    --max_cluster_bias_INS $MAX_CLUSTER_BIAS_INS \
    --diff_ratio_merging_INS $DIFF_RATIO_MERGING_INS \
    --max_cluster_bias_DEL $MAX_CLUSTER_BIAS_DEL \
    --diff_ratio_merging_DEL $DIFF_RATIO_MERGING_DEL \
    --genotype \
    --min_mapq 20 \
    ${input_bam} $ref38 ${outputname}.vcf $WORK_DIR --threads 8; } 2> "$TIME_LOG_TMP"
echo "cuteSV analysis complete. Main output: ${outputname}.vcf"


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
  echo "job_id: ${outputname}"
  echo "total_cpu_time: ${total_cpu_seconds}s (${user_seconds}s user + ${sys_seconds}s sys)"
  echo "wall_clock_time: ${real_time}"
  echo "---"
) >> running_time.log

# Clean up the temporary file.
rm "$TIME_LOG_TMP"
echo "Running time successfully logged to running_time.log"


# --- Post-processing ---
echo "Filtering for inversions..."
grep -E '^#|SVTYPE=INV' ${outputname}.vcf > inv_cutesv2_${outputname}.vcf
grep -v '#' inv_cutesv2_${outputname}.vcf > inv_woheader_cutesv2_${outputname}.vcf
echo "Inversion filtering complete."

