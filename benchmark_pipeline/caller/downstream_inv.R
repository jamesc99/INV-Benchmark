#!/usr/bin/env Rscript

# Load necessary libraries
suppressMessages({
  library(StructuralVariantAnnotation)
  library(VariantAnnotation)
  library(stringr)
})

# Command line arguments for the VCF file, genome, and output directory
args <- commandArgs(trailingOnly = TRUE)
vcf_file <- args[1]
genome <- args[2]
output_dir <- args[3]

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Read the VCF file
vcf <- readVcf(vcf_file, genome)

# Update header to include the SIMPLE_TYPE and SVLEN fields if not already present
info(header(vcf)) = unique(as(rbind(as.data.frame(info(header(vcf))), data.frame(
  row.names = c("SIMPLE_TYPE", "SVLEN"),
  Number = c("1", "1"),
  Type = c("String", "Integer"),
  Description = c("Simple event type annotation based purely on breakend position and orientation.",
                  "Length of the structural variant"))), "DataFrame"))

# Extract breakpoints and calculate SV type
breakpoints <- breakpointRanges(vcf)
svtype <- simpleEventType(breakpoints)

# Annotate the VCF with SIMPLE_TYPE and SVLEN
info(vcf)$SIMPLE_TYPE <- NA_character_
info(vcf)$SVLEN <- NA_integer_
info(vcf[breakpoints$sourceId])$SIMPLE_TYPE <- svtype
info(vcf[breakpoints$sourceId])$SVLEN <- breakpoints$svLen

# Write annotated VCF
annotated_vcf_file <- file.path(output_dir, "annotated.vcf")
writeVcf(vcf, annotated_vcf_file)

# Filter for high-confidence calls
high_confidence_gr <- breakpoints[breakpoints$FILTER == "PASS" & partner(breakpoints)$FILTER == "PASS"]

# Subset the VCF for high-confidence calls
high_confidence_vcf <- vcf[high_confidence_gr$sourceId]

# Write high-confidence annotated VCF
confident_vcf_file <- file.path(output_dir, "confident_annotated.vcf")
writeVcf(high_confidence_vcf, confident_vcf_file)

# Filter for inversion events from the high-confidence calls
inversion_gr <- high_confidence_gr[simpleEventType(high_confidence_gr) == "INV"]

# Create BED file for inversion events
simple_inv_bed <- data.frame(
  chrom = seqnames(inversion_gr),
  start = as.integer((start(inversion_gr) + end(inversion_gr)) / 2),
  end = as.integer((start(partner(inversion_gr)) + end(partner(inversion_gr))) / 2),
  name = simpleEventType(inversion_gr),
  score = inversion_gr$QUAL,
  strand = "."
)

# Remove duplicate breakends
simple_inv_bed <- simple_inv_bed[simple_inv_bed$start < simple_inv_bed$end, ]

# Write BED file containing only inversion events
simple_inv_bed_file <- file.path(output_dir, "simple_inv.bed")
write.table(simple_inv_bed, simple_inv_bed_file, quote = FALSE, sep = '\t', row.names = FALSE, col.names = FALSE)


