# Benchmark for simple and complex genome inversions

**This is the github repository to supplement our preprint : preprint link here**

The repository contains scripts and documentation for building and assessing the inversions benchmark set. It includes the script to build benchmark inversion set based on assembly files and Strand-seq signals, and the benchmarking process with seven SV callers on five samples (HG002, HG00733, HG02818, HG03486, NA19240).

There are two main components in this repository: Inversion benchmark build-up pipeline and Performance assessment.

## Data Availability
All VCF files used in this benchmark were uploaded to: zenodo link here.

## Assembly-based INV benchmark set build-up pipeline

**annotsv_individual.sh**

Description:
Runs AnnotSV
 on individual structural variant (SV) VCF files. Each SLURM array job processes one VCF file listed in vcflist.txt. It annotates SVs against the specified annotation directory and optionally produces both .tsv and .vcf outputs.

```
Inputs:

LIST_FILE — text file listing input VCF paths (default: vcflist.txt)

ANNOTSV_DIR — AnnotSV installation directory

ANN_DIR — AnnotSV annotation database directory

GENES_FILE (optional) — candidate gene list file

EMIT_VCF (flag) — whether to output annotated VCF (1) or not (0)

Outputs:

<OUT_DIR>/<SAMPLE>/
├─ <sample>.annotsv.tsv — annotated SV table
├─ <sample>.annotsv.vcf.gz (optional) — annotated VCF
```

