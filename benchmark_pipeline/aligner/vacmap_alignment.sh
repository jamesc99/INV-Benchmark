#!/bin/bash

#SBATCH --job-name=vacmap_align
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=4
#SBATCH --time=10-72:00:00
#SBATCH --partition=long
#SBATCH --mem=60G
#SBATCH -A proj-fs000f2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=siyuan.cheng@bcm.edu

#bamfile='/stornext/snfs4/next-gen/scratch/luis/hermann/data/bams/mosaic/hg0733/NHGRI_UCSC_panel-HG00733_63x.bam'
#bam_list="/users/u250191/ryan_scratch_ln/rotation_project/HG002_HG003_HG004/HG002_NA24385_son/minimap2/bam_name.list"

ulimit -n 4096 

ref37="/users/mmahmoud/home/public_workplace/scripts/snakefiles/test/hs37d5_mainchr.fa"
ref38="/stornext/snfs170/next-gen/scratch/ryan/reference/human-grch38.fasta"
ref38_decoy="/users/mmahmoud/home/source/hifi_panel/new_instalation/reference/GRCh38_masked_v2_decoy_gene.fasta" # from fixetflix
T2T="/users/sedlazec/mydir/projects/beta_globolin/demux_reads/analysis/reference/chm13v2.0.fasta"
RHESUS="/users/mmahmoud/home/projects/rhesus/reference/GCA_003339765.3_Mmul_10_genomic.maskPAR.withMT.mainchrs.fa"

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

samtools="/hgsc_software/samtools/samtools-1.19.2/bin/samtools"

. /users/u250191/.bashrc 

#mamba activate vacmap_env
mamba activate /hgsc_software/miniconda/miniconda3/envs/vacmap-1.0.1/

#vacmap -ref ${ref38} -read ${fastqfile} -mode S --MD -t 12 --rg-id ${rg_id} --rg-sm ${rg_sm} | ${samtools} sort -@12 > sorted_${fastq_basename}.bam
vacmap -ref ${ref38} -read ${fastqfile} -mode S --MD -t 4 --rg-id ${rg_id} --rg-sm ${rg_sm} > ${fastq_basename}.sam
#${samtools} sort -@12 ${fastq_basename}.sam > sorted_${fastq_basename}.bam
#vacmap -ref ${ref38} -read ${fastqfile} -mode S --MD -t 12 --rg-id ${rg_id} --rg-sm ${rg_sm} | ${samtools} sort -@12 > sorted_${fastq_basename}.bam
#${samtools} index -@12 sorted_${fastq_basename}.bam

#rm ${fastq_basename}.sam

#-ref The path of reference sequence.
#-read The path of long reads.
#-t The number of threads to use.
#-mode S For discovering complex variants (Pacbio CLR, ONT, HiFi).
#-mode H For aligning high error rate long read (Pacbio CLR, ONT).
#-mode L For aligning low error rate long read (Pacbio HiFi).
#--eqx Output =/X CIGAR operators for sequence match/mismatch.
#--MD Output the MD tag.
#--cs[=short|long] Output the cs tag. (deflaut: short cs).
#-k k-mer size (no larger than 28, deflaut: 15) # set -k 19 -w 10 for HiFi data  to reduce run time (2X faster) but there is very small decrease in accuracy.
#-w minimizer window size. (deflaut: 10)
#--rg-id <string>
#    Adds RG:Z:<string> to all alignments in SAM/BAM [none]
#--rg-sm <string>
#    RG header: Sample [none]


