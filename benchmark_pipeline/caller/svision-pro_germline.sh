#!/bin/bash
#SBATCH --job-name=svision-pro_germline
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=72:00:00
#SBATCH --partition=medium
#SBATCH -A proj-fs0002
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=siyuan.cheng@bcm.edu

# --- Environment Setup ---
source /users/u250191/.bashrc
conda activate svision-pro-env
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
ref38="/users/u250191/ryan_scratch_ln/reference/human-grch38.fasta"
access_bed="/users/u250191/ryan_scratch_ln/tools/SVision-pro-1.9/src/pre_process/hg38.access.10M.bed"
WORK_DIR=$(pwd)
default_model_path="/users/u250191/ryan_scratch_ln/tools/SVision-pro-1.9/src/pre_process/model_liteunet_256_8_16_32_32_32.pth"
bcftools=/hgsc_software/bcftools/bcftools-1.19/bin/bcftools

if [ -z "$1" ]; then
    echo "Error: No input BAM file provided."
    echo "Usage: sbatch $0 /path/to/your/file.bam"
    exit 1
fi
# WARNING: SVision-pro requires an absolute path for the input BAM file.
input_bam=$(realpath $1)


# --- Automatic Naming and Directory Setup ---
if [ -z "$SLURM_JOB_ID" ]; then
    echo "Error: This script must be run as a Slurm job." >&2
    exit 1
fi
job_id=$SLURM_JOB_ID

# Create a unique output directory for this job to prevent conflicts
OUT_DIR="${WORK_DIR}/svision_out_${job_id}"
mkdir -p "$OUT_DIR"
echo "Using Slurm Job ID ${job_id}. Output will be in: ${OUT_DIR}"


# --- Main Execution & Time Logging ---
TIME_LOG_TMP=$(mktemp)
echo "Starting SVision-pro analysis..."
# The time command's stderr is redirected to the temporary file.
{ time SVision-pro --out_path ${OUT_DIR} --target_path ${input_bam} \
	--genome_path ${ref38} --model_path ${default_model_path} --access_path ${access_bed} \
	--sample_name ${job_id} --detect_mode germline \
	--process_num 8; } 2> "$TIME_LOG_TMP"
echo "SVision-pro analysis complete."


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
echo "Running post-processing Python scripts..."
# Use the unique output directory for all file operations
python /users/u250191/ryan_scratch_ln/benchmark_inv/caller/svision-pro/extract_inv_add_simple_sv.py ${OUT_DIR}/${job_id}.vcf ${OUT_DIR}/inv_only_add_simplesv.vcf
grep -v '#' ${OUT_DIR}/inv_only_add_simplesv.vcf > ${OUT_DIR}/inv_only_add_simplesv_woheader.vcf

echo "Sorting VCF..."
$bcftools sort ${OUT_DIR}/inv_only_add_simplesv.vcf > ${OUT_DIR}/inv_only_add_simplesv.sorted.vcf
grep -v '#' ${OUT_DIR}/inv_only_add_simplesv.sorted.vcf > ${OUT_DIR}/inv_only_add_simplesv_woheader.sorted.vcf

# Clean up intermediate (unsorted) files
rm ${OUT_DIR}/inv_only_add_simplesv.vcf ${OUT_DIR}/inv_only_add_simplesv_woheader.vcf
echo "All steps complete."

