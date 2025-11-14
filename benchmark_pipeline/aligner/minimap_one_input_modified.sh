#!/bin/bash
#SBATCH --job-name=minimap2_mapping
#SBATCH --output=%x_%A_%a.out
#SBATCH --error=%x_%A_%a.err
#SBATCH --ntasks=6
#SBATCH --cpus-per-task=2
#SBATCH --time=10-72:00:00
#SBATCH --partition=long
#SBATCH -A proj-fs0002
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=siyuan.cheng@bcm.edu

# Reference genomes
ref37="/users/mmahmoud/home/public_workplace/scripts/snakefiles/test/hs37d5_mainchr.fa"
ref38="/users/u250191/ryan_scratch_ln/reference/human-grch38.fasta"
ref38_decoy="/users/mmahmoud/home/source/hifi_panel/new_instalation/reference/GRCh38_masked_v2_decoy_gene.fasta"
T2T="/users/sedlazec/mydir/projects/beta_globolin/demux_reads/analysis/reference/chm13v2.0.fasta"
RHESUS="/users/mmahmoud/home/projects/rhesus/reference/GCA_003339765.3_Mmul_10_genomic.maskPAR.withMT.mainchrs.fa"

minimap="/users/u250191/ryan_scratch_ln/tools/minimap2-2.28_x64-linux/minimap2"
samtools="/hgsc_software/samtools/samtools-1.19.2/bin/samtools"

# Get alignment type, FASTQ file path, and sample information from arguments
ALIGN_TYPE=$1
FASTQ_FILE=$2
INFOR=$3

# Check if ALIGN_TYPE, FASTQ_FILE, and INFOR are specified
if [ -z "$ALIGN_TYPE" ] || [ -z "$FASTQ_FILE" ] || [ -z "$INFOR" ]; then
  echo "Error: ALIGN_TYPE (ont or pacbio), FASTQ_FILE, and INFOR must be specified as arguments."
  exit 1
fi

# Extract the base name of the FASTQ file to use in output file names
BASE_NAME=$(basename "$FASTQ_FILE")

# Run minimap2 based on the specified alignment type and directly output BAM
case "$ALIGN_TYPE" in
  ont)
    ALIGN_OPTS="-ax map-ont"
    ;;
  pacbio)
    ALIGN_OPTS="-ax map-hifi -H"
    ;;
  *)
    echo "Error: Invalid ALIGN_TYPE specified. Use 'ont' or 'pacbio'."
    exit 1
    ;;
esac

INTERMEDIATE_BAM="intermediate_${BASE_NAME}.bam"
SORTED_BAM="sorted_${BASE_NAME}.bam"

$minimap -Y $ALIGN_OPTS -R "@RG\tSM:${INFOR}\tID:${INFOR}" "$ref38" "$FASTQ_FILE" --MD -t 12 > intermediate.sam
if [ $? -ne 0 ]; then
  echo "Error during minimap2 alignment for $ALIGN_TYPE."
  exit 1
fi

$samtools view -S -b -@ 12 intermediate.sam -o "$INTERMEDIATE_BAM"
if [ $? -ne 0 ]; then
  echo "Error converting SAM to BAM for $ALIGN_TYPE."
  exit 1
fi

$samtools sort -@ 12 -o "$SORTED_BAM" "$INTERMEDIATE_BAM"
if [ $? -ne 0 ]; then
  echo "Error sorting BAM for $ALIGN_TYPE."
  exit 1
fi

# Index the sorted BAM file
$samtools index -@ 12 "$SORTED_BAM"

echo "Pipeline completed successfully, and intermediate files have been retained."




