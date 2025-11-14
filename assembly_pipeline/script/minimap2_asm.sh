#!/bin/bash
#SBATCH --job-name=minimap2_asm
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

minimap2=/hgsc_software/minimap2/minimap2-2.22/bin/minimap2
samtools=/hgsc_software/samtools/samtools-1.21/bin/samtools
hg38_ref=/users/u250191/ryan_scratch_ln/reference/human-grch38.fasta

input_asm_fa=$1
temp_name=$2

$minimap2 -ax asm5 --eqx -t 16 $hg38_ref $input_asm_fa | $samtools view -@16 -b -o intermediate_${temp_name}.bam
$samtools sort -@16 -o sorted_${temp_name}.bam intermediate_${temp_name}.bam
$samtools index -@16 sorted_${temp_name}.bam

rm intermediate_${temp_name}.bam

