#!/bin/bash
#SBATCH --job-name=vacmap_asm
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=2
#SBATCH --mem=50gb
#SBATCH --time=72:00:00
#SBATCH --partition=medium
#SBATCH -A proj-fs0002
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=siyuan.cheng@bcm.edu

set -euo pipefail

# --- Binaries ---
VACMAP=vacmap   # or set absolute path if needed, e.g. /hgsc_software/vacmap/vacmap
samtools=/hgsc_software/samtools/samtools-1.21/bin/samtools

# --- Reference (gz as in your BAM header) ---
hg38_ref_gz=/users/u250191/ryan_scratch_ln/reference/human-grch38.fasta.gz

# --- Threads (8 tasks × 2 cpus-per-task = 16) ---
THREADS="${SLURM_CPUS_ON_NODE:-16}"

# --- Inputs ---
# 1) assembly FASTA (e.g., hg002_maternal.fasta)
# 2) temp/sample name    (e.g., hg002_maternal)
input_asm_fa=$1
temp_name=$2

workdir="./${temp_name}_temp"
mkdir -p "$workdir"

# --- Run VACmap → BAM, then sort & index ---
$VACMAP \
  -ref "$hg38_ref_gz" \
  -read "$input_asm_fa" \
  -mode asm \
  -t "$THREADS" \
  -workdir "$workdir" \
  --H --fakecigar --debug \
| $samtools view -@ "$THREADS" -b -o "intermediate_${temp_name}.bam" -

$samtools sort  -@ "$THREADS" -o "sorted_${temp_name}.bam" "intermediate_${temp_name}.bam"
$samtools index -@ "$THREADS" "sorted_${temp_name}.bam"

# --- Cleanup ---
rm -f "intermediate_${temp_name}.bam"
# keep workdir for debugging; remove if you want:
# rm -rf "$workdir"

