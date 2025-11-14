#!/bin/bash
#SBATCH --job-name=severus
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --mem=32gb
#SBATCH --time=72:00:00
#SBATCH --partition=medium

mamba activate severus_env

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

ref38="human-grch38.fasta"
tr_bed="/severus/human_GRCh38_no_alt_analysis_set.trf.bed"
WORK_DIR=$(pwd)
severus=/Severus-1.1/severus.py

if [ -z "$1" ]; then
    exit 1
fi
input_bam=$1

if [ -z "$SLURM_JOB_ID" ]; then
    exit 1
fi
job_id=$SLURM_JOB_ID

TIME_LOG_TMP=$(mktemp)

{ time $severus \
    --target-bam ${input_bam} \
    --out-dir ${WORK_DIR} \
    -t 12 \
    --vntr-bed ${tr_bed}; } 2> "$TIME_LOG_TMP"

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

grep -E '^#|SVTYPE=INV' all_SVs/severus_all.vcf > inv_severus.vcf
grep -v '#' inv_severus.vcf > inv_severus_woheader.vcf

mamba deactivate

python add_end_vcf_severus.py inv_severus.vcf
python severus_vcf2bed.py

