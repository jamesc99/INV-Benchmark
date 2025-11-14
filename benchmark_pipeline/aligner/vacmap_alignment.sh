#!/bin/bash

#SBATCH --job-name=vacmap_align
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=4
#SBATCH --time=10-72:00:00
#SBATCH --partition=long
#SBATCH --mem=60G

ulimit -n 4096 

ref38="/stornext/snfs170/next-gen/scratch/ryan/reference/human-grch38.fasta"

fastqfile=$1
rg_id=$2
rg_sm=$3

#--rg-id UNIQUE_ID --rg-sm SAMPLE_ID

if [[ "${fastqfile}" == *.gz ]]
then
    fastq_basename=$(basename "${fastqfile}" .fastq.gz)
else
    fastq_basename=$(basename "${fastqfile}" .fastq)
fi

samtools="/samtools-1.21/bin/samtools"

mamba activate vacmap-1.0.1/

vacmap -ref ${ref38} -read ${fastqfile} -mode S --MD -t 4 --rg-id ${rg_id} --rg-sm ${rg_sm} > ${fastq_basename}.sam
${samtools} sort -@12 ${fastq_basename}.sam > sorted_${fastq_basename}.bam
${samtools} index -@12 sorted_${fastq_basename}.bam

