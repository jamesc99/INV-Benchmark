# Benchmark for simple and complex genome inversions

**This is the github repository to supplement our preprint : preprint link here**

The repository contains scripts and documentation for building and assessing the inversions benchmark set. It includes the script to build benchmark inversion set, and the benchmarking process with seven SV callers on five samples (HG002, HG00733, HG02818, HG03486, NA19240).

There are two main components in this repository: inversion benchmark assembly pipeline and performance benchmark pipeline.

## Data Availability
All VCF files used in this benchmark were uploaded to: zenodo link here.

## Assembly-based INV set generation pipeline

**minimap2_asm.sh/vacmap_asm.sh**

Description:
Run minimap2/VACmap
 on assembly mode on haplotyped-resolved assembly FASTQ file to generate haplotype-specific BAM files.
```
Inputs:
INPUT_FASTA: Assembly FASTA file path (haplotype-specific)
OUTPUT_NAME: Prefix of output files


Outputs:
sorted_${OUTPUT_NAME}.bam: Assembly BAM file (haplotype-specific)
```

**asmBAM2vcf_integrated.py**

Description:
Main script to generate assembly-refined VCF based on haplotype-resolved assembly BAM files.
```
Inputs:
Srand-Seq BED: Strand-seq indicated inversion loci (provided in /assembly_pipeline/raw_strandseq_bed/)
Paternal_assembly_minimap2 BAM: Paternal assembly BAM file by minimap2
Maternal_assembly_minimap2 BAM: Maternal assembly BAM file by minimap2
Paternal_assembly_VACmap BAM: Paternal assembly BAM file by VACmap
Maternal_assembly_VACmap BAM: Maternal assembly BAM file by VACmap
TR_5Kb: Tandem Repeats annotation file for TR >= 5 Kb created by merging Adotto and UCSC Table Browser (provided in /assembly_pipeline/repeats_annotation/)
TR_10Kb: Tandem Repeats annotation file for TR >= 10 Kb created by merging Adotto and UCSC Table Browser (provided in /assembly_pipeline/repeats_annotation/)
Centromere: Centromere annotation file from UCSC Table Browser (provided in /assembly_pipeline/repeats_annotation/)
(optional) Debug: Debug mode by generating detailed log for each Strand-seq locus


Outputs:
Final VCF: Assembly-refined VCF file for each Strand-seq loci. Used for later benchmark (provided in /data_zenodo/truvari_ref_vcf/)
Dropped VCF: All inversion that were excluded from the final VCF, primarily those on chrY, in centromeres, or where both aligners reported NOSIGNAL
Statistics Report: A text report summarizing the contents of the Final VCF
Aligner Comparison Table: A comparison matrix showing how the two aligners agreed or disagreed before they were merged
```


## Benchmarking pipeline

### Step1: Aligner
**minimap2_alignment.sh**

Description:
Run minimap2 on alignment mode on FASTQ file to generate BAM files.
```
Inputs:
ALIGN_TYPE: Data type for sequencer-specific parameter ('pacbio' or 'ont')
FASTQ_FILE: Path for input fastq or fastq.gz file
INFOR: Information in RG tag

Outputs:
minimap2_BAM: Alignment BAM file
```

**vacmap_alignment.sh**

Description:
Run VACmap on alignment mode on FASTQ file to generate BAM files.
```
Inputs:
FASTQ_FILE: Path for input fastq or fastq.gz file
RG_ID: Information in RG_ID tag
RG_SM: Information in RG_SM tag

Outputs:
VACmap_BAM: Alignment BAM file
```

### Step2: Caller

**cutesv2_germline.sh**: Run cuteSV2 on alignment BAM file to generate INV calls.

**severus_default.sh**: Run Severus on alignment BAM file to generate INV calls.

**sniffles2.7_benchmark.sh**: Run Sniffles2 on alignment BAM file to generate INV calls.

**svim_alignment.sh**: Run SVIM on alignment BAM file to generate INV calls.

**svision-pro_germline.sh**: Run SVision-pro on alignment BAM file to generate INV calls.

**gridss2_germline_sr.sh**: Run GRIDSS2 on alignment BAM file to generate INV calls.

**manta_germline_shortread.sh**: Run Manta on alignment BAM file to generate INV calls.


### Step3: Truvari bench
**automated_truvari_bench.sh**

Description:
Run Truvari bench on caller-specific INV calls (VCF file) on INV trueset
```
Inputs:
INV_VCF: VCF file generated in Step2
Ref_name: Prefix used for this Truvari run

Outputs:
${Ref_name}_${base_name}: Folder contains Truvari bench output
```


































