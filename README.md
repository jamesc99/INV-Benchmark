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









