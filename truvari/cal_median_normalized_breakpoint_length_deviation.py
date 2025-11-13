#!/usr/bin/env python

import pysam
import sys
import statistics


def calculate_metrics(vcf_file):
    """
    Calculates absolute and normalized deviations for length and breakpoints
    for all INV variants in a VCF file.

    Args:
        vcf_file: An opened pysam.VariantFile object.

    Returns:
        A tuple containing four lists:
        - length_deviation (list of int): Absolute difference in length.
        - normalized_length_deviation (list of float): Length diff / true length.
        - breakpoint_deviation (list of int): Absolute difference in breakpoints.
        - normalized_breakpoint_deviation (list of float): Breakpoint diff / true length.
    """
    length_deviation = []
    normalized_length_deviation = []
    breakpoint_deviation = []
    normalized_breakpoint_deviation = []

    for record in vcf_file.fetch():
        if record.info.get('SVTYPE') == 'INV' and 'SVLEN' in record.info:
            # --- Get True Inversion Length ---
            svlen_val = record.info['SVLEN']
            if isinstance(svlen_val, (tuple, list)):
                inv_length = abs(svlen_val[0])
            else:
                inv_length = abs(svlen_val)

            if inv_length == 0:
                continue

            # --- Calculate Length Deviations ---
            if 'SizeDiff' in record.info:
                size_diff = abs(record.info['SizeDiff'])
                length_deviation.append(size_diff)
                normalized_length_deviation.append(size_diff / inv_length)

            # --- Calculate Breakpoint Deviations ---
            start_dev = abs(record.info.get('StartDistance', 0))
            end_dev = abs(record.info.get('EndDistance', 0))
            abs_bp_dev = start_dev + end_dev
            breakpoint_deviation.append(abs_bp_dev)
            normalized_breakpoint_deviation.append(abs_bp_dev / inv_length)

    return length_deviation, normalized_length_deviation, breakpoint_deviation, normalized_breakpoint_deviation


def main(vcf_path):
    """
    Main function to open VCF, run calculations, and write the report.
    """
    try:
        vcf_file = pysam.VariantFile(vcf_path)
    except (IOError, ValueError) as e:
        print(f"Error: Could not open or parse VCF file '{vcf_path}'.")
        print(f"Details: {e}")
        sys.exit(1)

    length_dev, norm_length_dev, abs_bp_dev, norm_bp_dev = calculate_metrics(vcf_file)

    if not length_dev:
        print("Warning: No valid INV variants with required INFO fields (SVTYPE, SVLEN, SizeDiff) found.")
        return

    output_log_file = "length_breakpoint_deviation.log"

    with open(output_log_file, "w") as log_file:
        log_file.write("Inversion Call Precision Analysis\n")
        log_file.write("=================================\n")
        log_file.write(f"Input VCF: {vcf_path}\n")
        log_file.write(f"Total INV variants processed: {len(length_dev)}\n\n")

        # Average and Median Absolute Length Deviation
        avg_length_dev = sum(length_dev) / len(length_dev)
        med_length_dev = statistics.median(length_dev)
        log_file.write(f"Average Absolute Length Deviation (SizeDiff): {avg_length_dev:.2f} bp\n")
        log_file.write(f"Median Absolute Length Deviation (SizeDiff): {med_length_dev:.2f} bp\n")

        # Average and Median Normalized Length Deviation
        avg_norm_length_dev = sum(norm_length_dev) / len(norm_length_dev)
        med_norm_length_dev = statistics.median(norm_length_dev)
        log_file.write(f"Average Normalized Length Deviation: {avg_norm_length_dev * 100:.4f}%\n")
        log_file.write(f"Median Normalized Length Deviation: {med_norm_length_dev * 100:.4f}%\n\n")

        # Average and Median Absolute Breakpoint Deviation
        avg_abs_bp_dev = sum(abs_bp_dev) / len(abs_bp_dev)
        med_abs_bp_dev = statistics.median(abs_bp_dev)
        log_file.write(f"Average Absolute Breakpoint Deviation: {avg_abs_bp_dev:.2f} bp\n")
        log_file.write(f"Median Absolute Breakpoint Deviation: {med_abs_bp_dev:.2f} bp\n")

        # Average and Median Normalized Breakpoint Deviation
        avg_norm_bp_dev = sum(norm_bp_dev) / len(norm_bp_dev)
        med_norm_bp_dev = statistics.median(norm_bp_dev)
        log_file.write(f"Average Normalized Breakpoint Deviation: {avg_norm_bp_dev * 100:.4f}%\n")
        log_file.write(f"Median Normalized Breakpoint Deviation: {med_norm_bp_dev * 100:.4f}%\n")

    print(f"Analysis complete. Results saved to {output_log_file}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python sv_precision_analysis.py <path/to/your/vcf_file.vcf.gz>")
        sys.exit(1)

    main(sys.argv[1])

