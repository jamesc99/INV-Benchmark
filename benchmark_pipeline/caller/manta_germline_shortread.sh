#!/bin/bash
#SBATCH --job-name=manta_germline
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --mem=32gb
#SBATCH --time=72:00:00
#SBATCH --partition=medium
#SBATCH -A proj-fs0002
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=siyuan.cheng@bcm.edu

source /users/u250191/.bashrc
mamba activate manta_env

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
    echo "$minutes * 60 + $seconds" | bc
}

MANTA_INSTALL_DIR="/hgsc_software/manta/manta-1.6.0"
REFERENCE_GENOME="/stornext/snfs170/next-gen/scratch/ryan/reference/human-grch38.fasta"
samtools="/hgsc_software/samtools/samtools-1.19.2/bin/samtools"
WORK_DIR=$(pwd)

if [ -z "$1" ]; then
    exit 1
fi
input_file=$1

if [ -z "$SLURM_JOB_ID" ]; then
    exit 1
fi
job_id=$SLURM_JOB_ID

# --- Manta Configuration & Time Logging (Step 1) ---
TIME_LOG_TMP1=$(mktemp)
{ time ${MANTA_INSTALL_DIR}/bin/configManta.py \
    --referenceFasta=${REFERENCE_GENOME} \
    --bam=${input_file} \
    --runDir=${WORK_DIR}; } 2> "$TIME_LOG_TMP1"

config_real_time_str=$(grep '^real' "$TIME_LOG_TMP1" | awk '{print $2}')
config_user_time_str=$(grep '^user' "$TIME_LOG_TMP1" | awk '{print $2}')
config_sys_time_str=$(grep '^sys' "$TIME_LOG_TMP1" | awk '{print $2}')
rm "$TIME_LOG_TMP1"

# --- Main Execution & Time Logging (Step 2) ---
TIME_LOG_TMP2=$(mktemp)
{ time ${WORK_DIR}/runWorkflow.py -j 8; } 2> "$TIME_LOG_TMP2"

workflow_real_time_str=$(grep '^real' "$TIME_LOG_TMP2" | awk '{print $2}')
workflow_user_time_str=$(grep '^user' "$TIME_LOG_TMP2" | awk '{print $2}')
workflow_sys_time_str=$(grep '^sys' "$TIME_LOG_TMP2" | awk '{print $2}')
rm "$TIME_LOG_TMP2"

# --- Process and Record Combined Running Time ---
config_real_seconds=$(time_to_seconds "$config_real_time_str")
config_user_seconds=$(time_to_seconds "$config_user_time_str")
config_sys_seconds=$(time_to_seconds "$config_sys_time_str")

workflow_real_seconds=$(time_to_seconds "$workflow_real_time_str")
workflow_user_seconds=$(time_to_seconds "$workflow_user_time_str")
workflow_sys_seconds=$(time_to_seconds "$workflow_sys_time_str")

total_wall_seconds=$(echo "$config_real_seconds + $workflow_real_seconds" | bc)
total_user_seconds=$(echo "$config_user_seconds + $workflow_user_seconds" | bc)
total_sys_seconds=$(echo "$config_sys_seconds + $workflow_sys_seconds" | bc)
total_cpu_seconds=$(echo "$total_user_seconds + $total_sys_seconds" | bc)

(
  echo "job_id: ${job_id}"
  echo "total_cpu_time: ${total_cpu_seconds}s (${total_user_seconds}s user + ${total_sys_seconds}s sys)"
  echo "wall_clock_time: ${total_wall_seconds}s"
  echo "---"
) >> running_time.log

# --- Post-processing ---
${MANTA_INSTALL_DIR}/libexec/convertInversion.py $samtools $REFERENCE_GENOME ./results/variants/diploidSV.vcf.gz > converted.vcf

grep -E '^#|SVTYPE=INV' converted.vcf > inv_with_infor.vcf
grep -v '^#' inv_with_infor.vcf > inv_only.vcf

awk 'BEGIN{OFS="\t"; print "Chr", "Start", "End", "Length", "GT"} {split($8, info, ";"); for(i in info) {split(info[i], kv, "="); if (kv[1] == "END") end_value = kv[2]; else if (kv[1] == "SVLEN") svlen_value = kv[2];} split($10, genotype_info, ":"); genotype = genotype_info[1]; print $1, $2, end_value, svlen_value, genotype;}' inv_only.vcf > inv_bedlike.bed

mamba deactivate

