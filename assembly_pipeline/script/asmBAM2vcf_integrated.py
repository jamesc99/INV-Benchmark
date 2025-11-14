#!/usr/bin/env python3

"""
Integrated Inversion Breakpoint Refinement Pipeline

This script combines 7 steps into a single executable pipeline to:
1.  Find candidate inversion breakpoints from paternal/maternal BAMs (for two aligners).
2.  Merge haplotype-specific VCFs into diploid calls.
3.  Filter out redundant/contained inversions based on the input BED.
4.  Annotate the VCFs with repeat and centromere information.
5.  Merge the results from the two aligners into a final, prioritized VCF.
6.  Generate statistics on the final merged VCF.
7.  Generate a comparison matrix for the two aligners' results before merging.
"""

# ############################################################################
# --- ALL IMPORTS ---
# ############################################################################

import sys
import os
import argparse
import statistics
import multiprocessing
import logging
import re
import collections
from collections import defaultdict

# Try to import required third-party libraries
try:
    import pysam
except ImportError:
    print("FATAL ERROR: `pysam` module not found. This pipeline requires `pysam`.", file=sys.stderr)
    print("Please install it: pip install pysam", file=sys.stderr)
    sys.exit(1)

try:
    import pandas as pd
except ImportError:
    print("FATAL ERROR: `pandas` module not found. This pipeline requires `pandas`.", file=sys.stderr)
    print("Please install it: pip install pandas", file=sys.stderr)
    sys.exit(1)

try:
    from intervaltree import IntervalTree
except ImportError:
    print("FATAL ERROR: `intervaltree` module not found. This pipeline requires `intervaltree`.", file=sys.stderr)
    print("Please install it: pip install intervaltree", file=sys.stderr)
    sys.exit(1)


# ############################################################################
# --- STEP 1: bam2vcfonbed_v5_multicore_debugLOG.py ---
# ############################################################################

# Get a logger instance for this module (used primarily by Step 1)
logger = logging.getLogger(__name__)

def setup_logging_step1(log_file="processing.log", debug_mode=False):
    """Sets up logging to file and console."""
    print(f"DEBUG_SETUP_LOGGING_ENTRY: Called with log_file='{log_file}', debug_mode={debug_mode}", flush=True)

    logger.setLevel(logging.DEBUG) # Set the logger to the lowest level, handlers will filter

    if logger.hasHandlers():
        logger.handlers.clear()

    # File handler
    try:
        fh = logging.FileHandler(log_file, mode='w')
        fh_formatter = logging.Formatter('%(asctime)s - %(levelname)s - [%(module)s.%(funcName)s:%(lineno)d] - %(message)s')
        fh.setFormatter(fh_formatter)
        fh.setLevel(logging.DEBUG if debug_mode else logging.INFO)
        logger.addHandler(fh)
        print(f"DEBUG_SETUP_LOGGING: File handler configured for '{log_file}'.", flush=True)
    except Exception as e:
        print(f"DEBUG_SETUP_LOGGING_ERROR: FAILED to configure file handler for '{log_file}': {e}", flush=True)
        with open("critical_file_handler_error.txt", "w") as f_crit_err:
            f_crit_err.write(f"Failed to configure file handler for '{log_file}': {e}\n")

    # Console handler
    try:
        ch = logging.StreamHandler(sys.stdout) # Log to standard output
        ch_formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        ch.setFormatter(ch_formatter)
        ch.setLevel(logging.INFO)
        if debug_mode:
            ch.setLevel(logging.DEBUG)
        logger.addHandler(ch)
        print(f"DEBUG_SETUP_LOGGING: Console handler configured.", flush=True)
    except Exception as e:
        print(f"DEBUG_SETUP_LOGGING_ERROR: FAILED to configure console handler: {e}", flush=True)
        with open("critical_console_handler_error.txt", "w") as f_crit_err:
            f_crit_err.write(f"Failed to configure console handler: {e}\n")

    if debug_mode:
        logger.info("Logger Test (INFO): DEBUG mode enabled. Detailed logs will be written to file and console.")
        logger.debug("Logger Test (DEBUG): This is a debug message after setup.")
    else:
        logger.info(f"Logger Test (INFO): Logging initialized. INFO logs to console. DEBUG/INFO to {log_file}.")
    logger.critical("Logger Test (CRITICAL): This is a critical message after setup. Should appear on console and file.")


def parse_bed_step1(bed_file):
    """
    Parse BED file and return list of inversion regions.
    (Step 1 Version)
    """
    inv_regions = []
    logger.info(f"Parsing BED file: {bed_file}")
    count = 0
    skipped_short = 0
    skipped_malformed = 0
    try:
        with open(bed_file, "r") as bed:
            for i, line in enumerate(bed):
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if len(fields) < 3:
                    logger.warning(f"Skipping malformed BED entry at line {i+1}: {line.strip()}")
                    skipped_malformed += 1
                    continue
                chrom = fields[0]
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError:
                    logger.warning(f"Skipping BED entry with non-integer coordinates at line {i+1}: {line.strip()}")
                    skipped_malformed += 1
                    continue

                inv_length = end - start
                if inv_length < 100:
                    skipped_short += 1
                    logger.debug(f"Skipping short inversion at line {i+1}: {chrom}:{start}-{end}, length {inv_length}bp")
                    continue

                window_extension = 3 * inv_length
                if window_extension > 3000000:
                    window_extension = 3000000

                search_window_start = max(0, start - window_extension)
                search_window_end = end + window_extension
                inv_regions.append((chrom, start, end, search_window_start, search_window_end, inv_length))
                logger.debug(f"[BED PARSE] Loaded: {chrom}:{start}-{end} (inv_length={inv_length}) "
                             f"=> search_window defined: {chrom}:{search_window_start}-{search_window_end}")
                count +=1
    except FileNotFoundError:
        logger.error(f"BED file not found: {bed_file}")
        print(f"FATAL_ERROR_PRINT: BED file not found: {bed_file}. Exiting.", flush=True)
        sys.exit(1)
    except Exception as e:
        logger.error(f"An error occurred while parsing BED file {bed_file}: {e}", exc_info=True)
        print(f"FATAL_ERROR_PRINT: Error parsing BED file {bed_file}: {e}. Exiting.", flush=True)
        sys.exit(1)

    logger.info(f"Parsed {count} inversion regions from BED file. Skipped {skipped_short} regions < 100bp. Skipped {skipped_malformed} malformed lines.")
    return inv_regions

def extract_candidate_positions_all_step1(read):
    """
    Extract candidate positions from a read using SA evidence or primary evidence.
    """
    candidates = []
    primary_strand = "+" if not read.is_reverse else "-"

    if read.mapq < 10: # Apply mapq filter upfront
        return candidates

    # Prioritize supplementary alignments (SA tag) if present
    if read.has_tag("SA"):
        sa_tag = read.get_tag("SA")
        for sa_alignment_str in sa_tag.split(";"):
            if not sa_alignment_str:
                continue
            fields = sa_alignment_str.split(",")
            if len(fields) < 6: # rname,pos,strand,CIGAR,mapQ,NM;
                logger.debug(f"Malformed SA entry for read {read.query_name}: {sa_alignment_str}")
                continue
            try:
                sa_pos = int(fields[1]) - 1  # 1-based in SA tag to 0-based
                sa_strand_char = fields[2]
                sa_mapq = int(fields[4])

                if sa_mapq >= 10:
                    candidates.append((sa_pos, sa_mapq, sa_strand_char))
            except ValueError as e:
                logger.debug(f"ValueError parsing SA tag field for read {read.query_name}: '{sa_alignment_str}'. Error: {e}")
                continue
    else:
        # Use primary alignment evidence if no SA tag is present
        if read.reference_start is not None: # Read is mapped
            candidates.append((read.reference_start, read.mapq, primary_strand))
        if read.reference_end is not None: # pysam's reference_end is exclusive
             candidates.append((read.reference_end -1, read.mapq, primary_strand)) # make it inclusive

        # Soft clipping
        if read.cigartuples:
            op_left, length_left = read.cigartuples[0]
            if op_left == 4 and length_left >= 10: # BAM CIGAR operation 'S'
                if read.reference_start is not None:
                    candidates.append((read.reference_start, read.mapq, primary_strand))

            op_right, length_right = read.cigartuples[-1]
            if op_right == 4 and length_right >= 10:
                 if read.reference_end is not None:
                    candidates.append((read.reference_end -1, read.mapq, primary_strand))

    return list(set(candidates))


def add_candidate_step1(candidate_dict, pos, mq, strand):
    if pos not in candidate_dict:
        candidate_dict[pos] = []
    candidate_dict[pos].append((mq, strand))

def candidate_search_all_step1(samfile, chrom, search_lower, search_upper):
    """
    Search for candidate positions in the window.
    """
    candidates = {}
    fetched_reads_count = 0
    processed_reads_count = 0
    logger.debug(f"Starting candidate search in {chrom}:{search_lower}-{search_upper}")
    try:
        for read_idx, read in enumerate(samfile.fetch(chrom, search_lower, search_upper)):
            fetched_reads_count += 1
            if read.is_unmapped or read.is_duplicate or read.is_qcfail or read.is_secondary:
                continue
            if read.mapq < 10:
                continue

            processed_reads_count +=1
            extracted_pos = extract_candidate_positions_all_step1(read)
            for pos, mq, strand in extracted_pos:
                if search_lower <= pos <= search_upper:
                    add_candidate_step1(candidates, pos, mq, strand)

            if fetched_reads_count > 0 and fetched_reads_count % 100000 == 0 :
                logger.info(f"[CANDIDATE SEARCH PROGRESS] {chrom}:{search_lower}-{search_upper} - Fetched {fetched_reads_count} reads, processed {processed_reads_count} valid reads. Current candidates: {len(candidates)}")

    except ValueError as e:
        logger.error(f"Error fetching reads for {chrom}:{search_lower}-{search_upper}. Error: {e}", exc_info=True)
        return {}
    except Exception as e:
        logger.error(f"Unexpected error during candidate search for {chrom}:{search_lower}-{search_upper}: {e}", exc_info=True)
        return {}

    logger.info(f"[CANDIDATE SEARCH DONE] {chrom}:{search_lower}-{search_upper} - Fetched {fetched_reads_count} reads, processed {processed_reads_count} valid reads. Found {len(candidates)} unique candidate positions.")
    return candidates

def cluster_candidates_step1(candidate_dict, cluster_threshold=100):
    """
    Cluster candidate positions.
    """
    if not candidate_dict:
        logger.debug("No candidates to cluster.")
        return {}

    positions = sorted(candidate_dict.keys())
    if not positions:
        return {}

    clusters = []
    current_cluster_positions = [positions[0]]
    for pos_idx in range(1, len(positions)):
        current_pos = positions[pos_idx]
        prev_pos_in_cluster = current_cluster_positions[-1]
        if current_pos - prev_pos_in_cluster <= cluster_threshold:
            current_cluster_positions.append(current_pos)
        else:
            clusters.append(list(current_cluster_positions))
            current_cluster_positions = [current_pos]

    if current_cluster_positions:
        clusters.append(list(current_cluster_positions))

    clustered_results = {}
    if not clusters:
        logger.debug("No clusters formed from candidates.")
        return {}

    for i, cluster_pos_list in enumerate(clusters):
        if not cluster_pos_list: continue

        try:
            consensus_pos = int(statistics.median(cluster_pos_list))
        except statistics.StatisticsError:
            logger.warning(f"Could not calculate median for cluster: {cluster_pos_list}. Skipping this cluster.")
            continue

        aggregated_support = []
        for pos_in_cluster in cluster_pos_list:
            aggregated_support.extend(candidate_dict[pos_in_cluster])

        if consensus_pos in clustered_results:
            logger.debug(f"Consensus position {consensus_pos} already exists. Merging support.")
            existing_support, existing_positions = clustered_results[consensus_pos]
            existing_support.extend(aggregated_support)
            existing_positions.extend(cluster_pos_list)
            clustered_results[consensus_pos] = (existing_support, sorted(list(set(existing_positions))))
        else:
            clustered_results[consensus_pos] = (aggregated_support, sorted(list(set(cluster_pos_list))))

        logger.debug(f"[CLUSTERING] Formed cluster {i+1}: original positions {cluster_pos_list} => consensus {consensus_pos} with {len(aggregated_support)} total supporting reads.")
    logger.info(f"Clustering resulted in {len(clustered_results)} final consensus clusters.")
    return clustered_results


def average_internal_distance_step1(positions_in_cluster):
    """
    Compute the average pairwise absolute difference.
    """
    unique_sorted_positions = sorted(list(set(positions_in_cluster)))
    n = len(unique_sorted_positions)
    if n < 2:
        logger.debug(f"Cluster {unique_sorted_positions} has < 2 unique positions, avg internal distance is 0.")
        return 0.0

    total_diff = 0
    pair_count = 0
    for i in range(n):
        for j in range(i + 1, n):
            total_diff += abs(unique_sorted_positions[j] - unique_sorted_positions[i])
            pair_count += 1

    avg_diff = total_diff / pair_count if pair_count > 0 else 0.0
    logger.debug(f"[AVG INTERNAL DISTANCE] Cluster footprint {unique_sorted_positions} (n={n}) => average distance {avg_diff:.2f}")
    return avg_diff

def select_candidate_pair_step1(start_clusters, end_clusters, bed_inv_length, bed_start, bed_end, chrom):
    """
    Evaluate candidate pairs from start and end clusters.
    """
    best_pair_info = None
    best_score = -float('inf')

    logger.debug(f"Starting candidate pair selection for inversion at {chrom}:{bed_start}-{bed_end} (length {bed_inv_length}bp).")
    logger.debug(f"Number of start clusters: {len(start_clusters)}, Number of end clusters: {len(end_clusters)}")

    if not start_clusters or not end_clusters:
        logger.info(f"Cannot select pair for {chrom}:{bed_start}-{bed_end}: No start or no end clusters found.")
        return None

    pair_eval_count = 0
    for s_consensus_pos, (s_support_list, s_cluster_positions) in start_clusters.items():
        s_total_support = len(s_support_list)
        if s_total_support == 0: continue

        s_plus_support = sum(1 for _, strand in s_support_list if strand == "+")
        s_minus_support = s_total_support - s_plus_support
        clarity_start = max(s_plus_support, s_minus_support) / s_total_support if s_total_support > 0 else 0
        dominant_start_strand = "+" if s_plus_support >= s_minus_support else "-"
        logger.debug(f"  Evaluating Start Cluster: consensus {s_consensus_pos}, support {s_total_support} ({s_plus_support}+, {s_minus_support}-), clarity {clarity_start:.2f}, dominant '{dominant_start_strand}'")

        for e_consensus_pos, (e_support_list, e_cluster_positions) in end_clusters.items():
            pair_eval_count += 1
            if e_consensus_pos <= s_consensus_pos:
                continue

            predicted_pair_length = e_consensus_pos - s_consensus_pos
            if predicted_pair_length < 100:
                logger.debug(f"    Skipping pair ({s_consensus_pos}, {e_consensus_pos}): length {predicted_pair_length}bp < 100bp.")
                continue

            e_total_support = len(e_support_list)
            if e_total_support == 0: continue

            e_plus_support = sum(1 for _, strand in e_support_list if strand == "+")
            e_minus_support = e_total_support - e_plus_support
            clarity_end = max(e_plus_support, e_minus_support) / e_total_support if e_total_support > 0 else 0
            dominant_end_strand = "+" if e_plus_support >= e_minus_support else "-"

            if dominant_start_strand == dominant_end_strand:
                logger.debug(f"    Skipping pair ({s_consensus_pos} ({dominant_start_strand}), {e_consensus_pos} ({dominant_end_strand})): Dominant strands are the same.")
                continue

            logger.debug(f"    Considering End Cluster: consensus {e_consensus_pos}, support {e_total_support} ({e_plus_support}+, {e_minus_support}-), clarity {clarity_end:.2f}, dominant '{dominant_end_strand}'")

            overall_clarity_score = min(clarity_start, clarity_end)
            deviation_from_bed = abs(s_consensus_pos - bed_start) + abs(e_consensus_pos - bed_end)
            norm_factor_prox = 4 * bed_inv_length if bed_inv_length > 0 else 1000000
            proximity_score = 1.0 - (deviation_from_bed / norm_factor_prox)
            proximity_score = max(0, min(proximity_score, 1))

            total_pair_support = s_total_support + e_total_support
            support_score = min(total_pair_support / 60.0, 1.0)

            length_similarity_score = 0.0
            if bed_inv_length > 0:
                if predicted_pair_length >= bed_inv_length:
                    diff = predicted_pair_length - bed_inv_length
                    if diff <= bed_inv_length:
                        length_similarity_score = 1.0
                    elif diff <= 3 * bed_inv_length:
                        length_similarity_score = 1.0 - ((diff - bed_inv_length) / (2 * bed_inv_length))
                    else:
                        length_similarity_score = 0.0
                else: # predicted_pair_length < bed_inv_length
                    if predicted_pair_length >= 0.5 * bed_inv_length:
                        length_similarity_score = 1.0
                    elif predicted_pair_length >= 0.25 * bed_inv_length:
                        length_similarity_score = 1.0 - ((0.5 * bed_inv_length - predicted_pair_length) / (0.25 * bed_inv_length))
                    else:
                        length_similarity_score = 0.0
            elif predicted_pair_length == 0 :
                 length_similarity_score = 1.0
            length_similarity_score = max(0, min(length_similarity_score, 1))


            current_composite_score = (
                0.60 * proximity_score +
                0.10 * overall_clarity_score +
                0.15 * support_score +
                0.15 * length_similarity_score
            )

            logger.debug(f"      Pair ({s_consensus_pos}, {e_consensus_pos}): Length {predicted_pair_length}bp. BED Length: {bed_inv_length}bp.")
            logger.debug(f"        Scores: Proximity={proximity_score:.3f} (deviation={deviation_from_bed}), Clarity={overall_clarity_score:.3f}, Support={support_score:.3f} (total reads={total_pair_support}), LengthSim={length_similarity_score:.3f}")
            logger.debug(f"        Composite Score: {current_composite_score:.4f}")

            if current_composite_score > best_score:
                best_score = current_composite_score
                best_pair_info = (s_consensus_pos, s_support_list, s_cluster_positions,
                                  e_consensus_pos, e_support_list, e_cluster_positions,
                                  current_composite_score, predicted_pair_length)
                logger.debug(f"        ----> NEW BEST PAIR for {chrom}:{bed_start}-{bed_end} with score {best_score:.4f}")

    if best_pair_info:
        logger.info(f"Selected best pair for {chrom}:{bed_start}-{bed_end} -> Start:{best_pair_info[0]}, End:{best_pair_info[3]}, Score:{best_pair_info[6]:.4f}, Length:{best_pair_info[7]}bp")
    else:
        logger.info(f"No suitable candidate pair found for {chrom}:{bed_start}-{bed_end} after evaluating {pair_eval_count} potential pairs.")

    return best_pair_info


def find_breakpoints_step1(bam_file_path, list_of_inv_regions, source_allele):
    """
    Main logic for finding breakpoints for a list of inversion regions.
    """
    samfile = None
    try:
        logger.debug(f"Opening BAM file: {bam_file_path}")
        samfile = pysam.AlignmentFile(bam_file_path, "rb")
        if not samfile.header:
             logger.error(f"BAM file {bam_file_path} has no header or is invalid. Cannot process.")
             print(f"FATAL_ERROR_PRINT: BAM file {bam_file_path} has no header or is invalid. Cannot process. For source {source_allele}", flush=True)
             return []
    except Exception as e:
        logger.error(f"Could not open or validate BAM file {bam_file_path}: {e}", exc_info=True)
        print(f"FATAL_ERROR_PRINT: Could not open or validate BAM file {bam_file_path}: {e}. For source {source_allele}", flush=True)
        return []

    vcf_output_entries = []
    num_total_regions = len(list_of_inv_regions)
    logger.info(f"Processing {num_total_regions} regions for source '{source_allele}' from BAM '{bam_file_path}'")

    for region_idx, inv_region_details in enumerate(list_of_inv_regions):
        chrom, bed_start_coord, bed_end_coord, global_search_start, global_search_end, bed_inv_len = inv_region_details

        logger.info(f"--- Processing Region {region_idx+1}/{num_total_regions} ({source_allele}): {chrom}:{bed_start_coord}-{bed_end_coord} (len: {bed_inv_len}bp) ---")

        flank_search_span = 3 * bed_inv_len
        if flank_search_span > 3000000: flank_search_span = 3000000

        start_bp_search_window_lower = max(0, bed_start_coord - flank_search_span)
        start_bp_search_window_upper = bed_start_coord + flank_search_span

        end_bp_search_window_lower = max(0, bed_end_coord - flank_search_span)
        end_bp_search_window_upper = bed_end_coord + flank_search_span

        logger.debug(f"Start breakpoint candidate search window: {chrom}:{start_bp_search_window_lower}-{start_bp_search_window_upper}")
        start_candidate_positions = candidate_search_all_step1(samfile, chrom, start_bp_search_window_lower, start_bp_search_window_upper)

        logger.debug(f"End breakpoint candidate search window: {chrom}:{end_bp_search_window_lower}-{end_bp_search_window_upper}")
        end_candidate_positions = candidate_search_all_step1(samfile, chrom, end_bp_search_window_lower, end_bp_search_window_upper)

        start_clusters_dict = cluster_candidates_step1(start_candidate_positions, cluster_threshold=100)
        end_clusters_dict = cluster_candidates_step1(end_candidate_positions, cluster_threshold=100)
        logger.info(f"Found {len(start_clusters_dict)} start clusters and {len(end_clusters_dict)} end clusters for {chrom}:{bed_start_coord}-{bed_end_coord}.")

        if not start_clusters_dict or not end_clusters_dict:
            logger.warning(f"No candidate clusters on one or both sides for {chrom}:{bed_start_coord}-{bed_end_coord}. Labeling as NOSIGNAL.")
            info_str = f"SVTYPE=INV;END={bed_end_coord + 1};LABEL=no_signal;SOURCE={source_allele};ORIGINAL_LOCUS={chrom}_{bed_start_coord}_{bed_end_coord};SVLEN={bed_inv_len}"
            vcf_output_entries.append(f"{chrom}\t{bed_start_coord + 1}\t.\tN\t<INV>\t.\tNOSIGNAL\t{info_str}")
            continue

        best_selected_pair = select_candidate_pair_step1(start_clusters_dict, end_clusters_dict, bed_inv_len, bed_start_coord, bed_end_coord, chrom)

        if best_selected_pair is None:
            logger.warning(f"No suitable candidate pair found for {chrom}:{bed_start_coord}-{bed_end_coord} after scoring. Labeling as NOSIGNAL.")
            info_str = f"SVTYPE=INV;END={bed_end_coord + 1};LABEL=no_signal;SOURCE={source_allele};ORIGINAL_LOCUS={chrom}_{bed_start_coord}_{bed_end_coord};SVLEN={bed_inv_len}"
            vcf_output_entries.append(f"{chrom}\t{bed_start_coord + 1}\t.\tN\t<INV>\t.\tNOSIGNAL\t{info_str}")
            continue

        final_s_pos, s_sup_reads, s_raw_positions_list, final_e_pos, e_sup_reads, e_raw_positions_list, pair_score, refined_svlen = best_selected_pair
        
        s_raw_positions_sorted_unique = sorted(list(set(s_raw_positions_list)))
        e_raw_positions_sorted_unique = sorted(list(set(e_raw_positions_list)))

        cipos_val1 = min(s_raw_positions_sorted_unique) - bed_start_coord if s_raw_positions_sorted_unique else 0
        cipos_val2 = max(s_raw_positions_sorted_unique) - bed_start_coord if s_raw_positions_sorted_unique else 0
        cipos_str_v1 = f"{cipos_val1},{cipos_val2}"

        ciend_val1 = min(e_raw_positions_sorted_unique) - bed_end_coord if e_raw_positions_sorted_unique else 0
        ciend_val2 = max(e_raw_positions_sorted_unique) - bed_end_coord if e_raw_positions_sorted_unique else 0
        ciend_str_v1 = f"{ciend_val1},{ciend_val2}"

        avg_id_start = average_internal_distance_step1(s_raw_positions_list) if len(set(s_raw_positions_list)) >= 5 else 0.0
        avg_id_end = average_internal_distance_step1(e_raw_positions_list) if len(set(e_raw_positions_list)) >= 5 else 0.0
        
        num_valid_opp_strand_pairs = 0
        unique_s_clusters_in_valid_pairs = set()
        unique_e_clusters_in_valid_pairs = set()

        for s_c, (s_sup, _) in start_clusters_dict.items():
            s_plus = sum(1 for _, s_str_val in s_sup if s_str_val == '+')
            dom_s = '+' if s_plus >= (len(s_sup) - s_plus) else '-'
            for e_c, (e_sup, _) in end_clusters_dict.items():
                if e_c <= s_c: continue
                e_plus = sum(1 for _, e_str_val in e_sup if e_str_val == '+')
                dom_e = '+' if e_plus >= (len(e_sup) - e_plus) else '-'
                if dom_s != dom_e:
                    num_valid_opp_strand_pairs += 1
                    unique_s_clusters_in_valid_pairs.add(s_c)
                    unique_e_clusters_in_valid_pairs.add(e_c)
        
        len_valid_uniq_s_clusters_v1 = len(unique_s_clusters_in_valid_pairs)
        len_valid_uniq_e_clusters_v1 = len(unique_e_clusters_in_valid_pairs)
        valid_pair_count_v1 = num_valid_opp_strand_pairs

        range_s_cluster = max(s_raw_positions_list) - min(s_raw_positions_list) if s_raw_positions_list else 0
        range_e_cluster = max(e_raw_positions_list) - min(e_raw_positions_list) if e_raw_positions_list else 0
        precision_ratio = 0.0
        if refined_svlen > 0:
            precision_ratio = max(range_s_cluster, range_e_cluster) / float(refined_svlen)
        elif max(range_s_cluster, range_e_cluster) > 0 :
            precision_ratio = 1.0
        
        filter_tag = "LOWCONF"
        filter_reason = "Default, or ratio > 0.50"

        if (len(unique_s_clusters_in_valid_pairs) >= 4 or len(unique_e_clusters_in_valid_pairs) >= 4) or \
           (num_valid_opp_strand_pairs >= 5):
            filter_tag = "PASS_IMPRECISE"
            filter_reason = "High ambiguity: >=4 unique clusters on one side or >=5 valid pairs"
        elif (avg_id_start > 10.0 and len(set(s_raw_positions_list)) >= 5) or \
             (avg_id_end > 10.0 and len(set(e_raw_positions_list)) >= 5):
            filter_tag = "PASS_IMPRECISE"
            filter_reason = f"High internal distance in chosen cluster (start_avg_id:{avg_id_start:.2f}, end_avg_id:{avg_id_end:.2f})"
        elif precision_ratio <= 0.20:
            filter_tag = "PASS_PRECISE"
            filter_reason = f"Precise: range/svlen ratio ({precision_ratio:.3f}) <= 0.20"
        elif precision_ratio <= 0.50:
            filter_tag = "PASS_IMPRECISE"
            filter_reason = f"Imprecise: range/svlen ratio ({precision_ratio:.3f}) > 0.20 and <= 0.50"
        
        logger.info(f"  FILTER decision for {chrom}:{final_s_pos}-{final_e_pos} ({source_allele}): {filter_tag} (Reason: {filter_reason})")

        label_v1 = "INV"
        final_end_v1 = final_e_pos + 1

        info_list_v1 = [
            f"SVTYPE=INV",
            f"END={final_end_v1}",
            f"LABEL={label_v1}",
            f"CIPOS={cipos_str_v1}",
            f"CIEND={ciend_str_v1}",
            f"SOURCE={source_allele}",
            f"ORIGINAL_LOCUS={chrom}_{bed_start_coord}_{bed_end_coord}",
            f"SVLEN={refined_svlen}",
            f"SUPMQ_START={','.join(str(mq) for mq, _ in s_sup_reads)}",
            f"SUPSTRAND_START={','.join(s_str for _, s_str in s_sup_reads)}",
            f"SUPMQ_END={','.join(str(mq) for mq, _ in e_sup_reads)}",
            f"SUPSTRAND_END={','.join(e_str for _, e_str in e_sup_reads)}",
            f"VALID_UNIQ_CLUSTER_START={len_valid_uniq_s_clusters_v1}",
            f"VALID_UNIQ_CLUSTER_END={len_valid_uniq_e_clusters_v1}",
            f"VALID_CLUSTER_PAIR={valid_pair_count_v1}"
        ]
        info_str = ";".join(info_list_v1)

        vcf_line = f"{chrom}\t{final_s_pos + 1}\t.\tN\t<INV>\t.\t{filter_tag}\t{info_str}" # VCF POS is 1-based
        vcf_output_entries.append(vcf_line)
        logger.debug(f"Added VCF entry: {vcf_line}")

    if samfile:
        logger.debug(f"Closing BAM file: {bam_file_path}")
        samfile.close()

    logger.info(f"Finished processing all regions for source '{source_allele}' from BAM '{bam_file_path}'. Generated {len(vcf_output_entries)} VCF entries.")
    return vcf_output_entries


def find_breakpoints_for_region_mp_wrapper_step1(args_bundle):
    """
    Wrapper for multiprocessing: unpacks arguments and calls find_breakpoints.
    """
    bam_file_path, single_region_data_list, source_allele, log_file_path, is_debug_mode = args_bundle
    process_id = multiprocessing.current_process().name

    logger.info(f"[{process_id}] Worker started for region: {single_region_data_list[0][0]}:{single_region_data_list[0][1]}-{single_region_data_list[0][2]} ({source_allele})")
    results = find_breakpoints_step1(bam_file_path, single_region_data_list, source_allele)
    logger.info(f"[{process_id}] Worker finished for region. Found {len(results)} entries.")
    return results


def _write_vcf_headers_step1(vcf_file_handle):
    """Helper function to write VCF headers, matching v1_old format."""
    logger.debug("Writing VCF headers (v1_old style).")
    vcf_file_handle.write("##fileformat=VCFv4.2\n")
    vcf_file_handle.write("##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Structural variant type\">\n")
    vcf_file_handle.write("##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of structural variant\">\n")
    vcf_file_handle.write("##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Variant length calculated from breakpoints\">\n")
    vcf_file_handle.write("##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description=\"Imprecise variant\">\n")
    vcf_file_handle.write("##INFO=<ID=CIPOS,Number=1,Type=String,Description=\"Confidence interval around POS for imprecise inversions, relative to original BED start (min_offset,max_offset)\">\n")
    vcf_file_handle.write("##INFO=<ID=CIEND,Number=1,Type=String,Description=\"Confidence interval around END for imprecise inversions, relative to original BED end (min_offset,max_offset)\">\n")
    vcf_file_handle.write("##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Indicates if variant is from paternal or maternal BAM\">\n")
    vcf_file_handle.write("##INFO=<ID=LABEL,Number=1,Type=String,Description=\"Structural variant label\">\n")
    vcf_file_handle.write("##INFO=<ID=ORIGINAL_LOCUS,Number=1,Type=String,Description=\"Original inversion coordinates from BED file\">\n")
    vcf_file_handle.write("##INFO=<ID=SUPMQ_START,Number=.,Type=String,Description=\"Mapping qualities of supporting reads at the START boundary\">\n")
    vcf_file_handle.write("##INFO=<ID=SUPSTRAND_START,Number=.,Type=String,Description=\"Strand information of supporting reads at the START boundary\">\n")
    vcf_file_handle.write("##INFO=<ID=SUPMQ_END,Number=.,Type=String,Description=\"Mapping qualities of supporting reads at the END boundary\">\n")
    vcf_file_handle.write("##INFO=<ID=SUPSTRAND_END,Number=.,Type=String,Description=\"Strand information of supporting reads at the END boundary\">\n")
    vcf_file_handle.write("##INFO=<ID=VALID_UNIQ_CLUSTER_START,Number=1,Type=Integer,Description=\"Number of unique start clusters involved in valid pairs (v1 style)\">\n")
    vcf_file_handle.write("##INFO=<ID=VALID_UNIQ_CLUSTER_END,Number=1,Type=Integer,Description=\"Number of unique end clusters involved in valid pairs (v1 style)\">\n")
    vcf_file_handle.write("##INFO=<ID=VALID_CLUSTER_PAIR,Number=1,Type=Integer,Description=\"Total number of potential valid (opposite strand) cluster pairs (v1 style)\">\n")
    vcf_file_handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def write_vcf_step1(output_vcf_path, paternal_entries_list, maternal_entries_list):
    logger.info(f"Starting VCF writing process for: {output_vcf_path}")

    def get_chrom_sort_key(chrom_name_str):
        c = chrom_name_str.lower().replace('chr', '')
        if c.isdigit(): return int(c)
        if c == 'x': return 23
        if c == 'y': return 24
        if c == 'm' or c == 'mt': return 25
        return 26

    parsed_vcf_entries = []
    all_raw_entries = paternal_entries_list + maternal_entries_list
    logger.debug(f"Total raw entries to process for VCF: {len(all_raw_entries)}")

    for vcf_line_str in all_raw_entries:
        fields = vcf_line_str.strip().split('\t')
        if len(fields) < 8:
            logger.warning(f"Skipping malformed VCF line during parsing: {vcf_line_str}")
            continue

        chrom, pos_str, _, _, _, _, fltr_str, info_str = fields[:8]
        pos = int(pos_str)

        original_locus_val = ""
        source_val = ""
        info_parts = info_str.split(';')
        for item in info_parts:
            if item.startswith("ORIGINAL_LOCUS="):
                original_locus_val = item.split("=",1)[1]
            elif item.startswith("SOURCE="):
                source_val = item.split("=",1)[1]
        
        if not original_locus_val:
            logger.warning(f"Entry {chrom}:{pos} is missing ORIGINAL_LOCUS in INFO. Using placeholder.")
            original_locus_val = f"MISSINGLOCUS_{chrom}_{pos}_{source_val}"

        parsed_vcf_entries.append({
            'line': vcf_line_str, 'chrom': chrom, 'pos': pos,
            'filter': fltr_str, 'info_str': info_str,
            'original_locus': original_locus_val, 'source': source_val
        })
    
    grouped_by_locus = defaultdict(list)
    for entry_dict in parsed_vcf_entries:
        grouped_by_locus[entry_dict['original_locus']].append(entry_dict)

    logger.info(f"Homogenizing filters across {len(grouped_by_locus)} original loci.")
    for locus, entries_for_locus in grouped_by_locus.items():
        if locus.startswith("MISSINGLOCUS_"): continue

        if len(entries_for_locus) == 2:
            entry1 = entries_for_locus[0]
            entry2 = entries_for_locus[1]
            
            sources_in_pair = {entry1.get('source', '').upper(), entry2.get('source', '').upper()}
            if not ({'PATERNAL', 'MATERNAL'} == sources_in_pair):
                logger.debug(f"Locus {locus} has two entries but not one PATERNAL and one MATERNAL. Sources: {sources_in_pair}. Skipping filter homogenization.")
                continue

            filter1 = entry1['filter']
            filter2 = entry2['filter']

            if ('PASS_PRECISE' in filter1 and 'PASS_IMPRECISE' in filter2) or \
               ('PASS_IMPRECISE' in filter1 and 'PASS_PRECISE' in filter2):
                logger.debug(f"Homogenizing filter for locus {locus}: {filter1}/{filter2} -> PASS_IMPRECISE for both.")
                
                line_parts1 = entry1['line'].split('\t')
                line_parts1[6] = 'PASS_IMPRECISE'
                entry1['line'] = '\t'.join(line_parts1)
                entry1['filter'] = 'PASS_IMPRECISE'

                line_parts2 = entry2['line'].split('\t')
                line_parts2[6] = 'PASS_IMPRECISE'
                entry2['line'] = '\t'.join(line_parts2)
                entry2['filter'] = 'PASS_IMPRECISE'

    final_entries_to_sort = []
    for locus, entries_for_locus in grouped_by_locus.items():
        final_entries_to_sort.extend(entries_for_locus)

    def get_sort_key_for_vcf_entry_v1_style(entry_dict):
        chrom_key = get_chrom_sort_key(entry_dict['chrom'])
        bed_s, bed_e = 0, 0
        
        original_locus_str = entry_dict.get('original_locus', '')
        if original_locus_str and not original_locus_str.startswith("MISSINGLOCUS_"):
            try:
                parts = original_locus_str.split('_')
                if len(parts) >= 3:
                    bed_s = int(parts[-2])
                    bed_e = int(parts[-1])
            except ValueError:
                logger.warning(f"Could not parse bed_start/end from ORIGINAL_LOCUS '{original_locus_str}' for sorting. Entry: {entry_dict['line']}")
            except IndexError:
                 logger.warning(f"Index error parsing bed_start/end from ORIGINAL_LOCUS '{original_locus_str}' for sorting. Entry: {entry_dict['line']}")

        source_priority = 0 if entry_dict.get('source', '').upper() == 'MATERNAL' else 1
        return (chrom_key, bed_s, bed_e, source_priority, entry_dict['pos'])

    logger.info(f"Sorting {len(final_entries_to_sort)} VCF entries using v1-style key.")
    try:
        sorted_final_entries = sorted(final_entries_to_sort, key=get_sort_key_for_vcf_entry_v1_style)
    except Exception as e:
        logger.error(f"Error during VCF entry sorting: {e}. Writing unsorted entries if any.", exc_info=True)
        sorted_final_entries = final_entries_to_sort

    try:
        with open(output_vcf_path, "w") as vcf_out_f:
            _write_vcf_headers_step1(vcf_out_f)
            for entry_dict in sorted_final_entries:
                vcf_out_f.write(entry_dict['line'] + "\n")
        logger.info(f"Successfully wrote {len(sorted_final_entries)} entries to VCF file: {output_vcf_path}")
    except IOError as e:
        logger.error(f"Could not write to VCF file {output_vcf_path}: {e}", exc_info=True)


def run_step1(args_dict):
    """
    Main function for Step 1, adapted to run from an arguments dictionary.
    """
    # Simulate the args object
    args = argparse.Namespace(**args_dict)

    try:
        setup_logging_step1(log_file=args.log_file, debug_mode=args.debug)
    except Exception as e:
        print(f"DEBUG_MAIN_LOGGING_SETUP_FAIL: CRITICAL - setup_logging failed: {e}", flush=True)
        with open("setup_logging_failed_at_step1.txt", "w") as f_err:
            f_err.write(f"setup_logging failed: {e}\nLog file was to be: {args.log_file}\nDebug mode: {args.debug}\n")
        sys.exit(1)

    logger.info("--- Inversion Breakpoint Finder Script Started (Official Log Start) ---")
    logger.info(f"Run parameters: BED='{args.bed}', PaternalBAM='{args.paternal_bam}', MaternalBAM='{args.maternal_bam}', OutputVCF='{args.output_vcf}', Threads={args.threads}, DebugMode={args.debug}, LogFile='{args.log_file}'")

    all_inv_regions = parse_bed_step1(args.bed)
    if not all_inv_regions:
        logger.error("No inversion regions were successfully parsed from the BED file. Exiting.")
        return

    num_cpus_to_use = args.threads
    logger.info(f"Will process {len(all_inv_regions)} regions using up to {num_cpus_to_use} thread(s).")

    paternal_vcf_results = []
    maternal_vcf_results = []

    if num_cpus_to_use > 1 and len(all_inv_regions) > 0:
        logger.info("Initializing multiprocessing pool.")
        paternal_mp_args = [(args.paternal_bam, [region_data], "PATERNAL", args.log_file, args.debug) for region_data in all_inv_regions]
        maternal_mp_args = [(args.maternal_bam, [region_data], "MATERNAL", args.log_file, args.debug) for region_data in all_inv_regions]
        
        paternal_results_ok = False
        maternal_results_ok = False
        try:
            with multiprocessing.Pool(processes=num_cpus_to_use) as pool:
                logger.info(f"Processing paternal BAM regions with {num_cpus_to_use} workers...")
                paternal_mp_outputs = pool.map(find_breakpoints_for_region_mp_wrapper_step1, paternal_mp_args)
                paternal_vcf_results = [entry for sublist_result in paternal_mp_outputs for entry in sublist_result if sublist_result]
                paternal_results_ok = True
                logger.info(f"Completed paternal BAM processing. Got {len(paternal_vcf_results)} entries.")

                logger.info(f"Processing maternal BAM regions with {num_cpus_to_use} workers...")
                maternal_mp_outputs = pool.map(find_breakpoints_for_region_mp_wrapper_step1, maternal_mp_args)
                maternal_vcf_results = [entry for sublist_result in maternal_mp_outputs for entry in sublist_result if sublist_result]
                maternal_results_ok = True
                logger.info(f"Completed maternal BAM processing. Got {len(maternal_vcf_results)} entries.")
        except Exception as e:
            logger.error(f"An error occurred during multiprocessing: {e}", exc_info=True)
            if not paternal_results_ok: paternal_vcf_results = []
            if not maternal_results_ok: maternal_vcf_results = []
            
            if not (paternal_vcf_results or maternal_vcf_results):
                 logger.error("Multiprocessing failed to produce any results. Will attempt sequential fallback.")
                 num_cpus_to_use = 1
            else:
                 logger.warning("Multiprocessing may have partially failed. Proceeding with any obtained results.")
    
    run_sequentially = False
    if num_cpus_to_use <= 1:
        run_sequentially = True
    elif not all_inv_regions:
        logger.info("No regions to process.")
    elif not (paternal_vcf_results or maternal_vcf_results) and args.threads > 1:
        logger.info("Falling back to sequential processing as multiprocessing yielded no results.")
        run_sequentially = True

    if run_sequentially and len(all_inv_regions) > 0:
        if args.threads > 1 and not (paternal_vcf_results or maternal_vcf_results) :
             pass
        else:
            logger.info("Starting sequential processing.")

        if not paternal_vcf_results:
            logger.info("Processing paternal BAM regions sequentially...")
            paternal_vcf_results = find_breakpoints_step1(args.paternal_bam, all_inv_regions, "PATERNAL")
            logger.info(f"Completed paternal BAM processing sequentially. Got {len(paternal_vcf_results)} entries.")
        else:
            logger.info("Skipping sequential paternal processing, results already exist.")

        if not maternal_vcf_results:
            logger.info("Processing maternal BAM regions sequentially...")
            maternal_vcf_results = find_breakpoints_step1(args.maternal_bam, all_inv_regions, "MATERNAL")
            logger.info(f"Completed maternal BAM processing sequentially. Got {len(maternal_vcf_results)} entries.")
        else:
            logger.info("Skipping sequential maternal processing, results already exist.")

    write_vcf_step1(args.output_vcf, paternal_vcf_results, maternal_vcf_results)
    logger.info("--- Inversion Breakpoint Finder Script (Step 1) Finished ---")


# ############################################################################
# --- STEP 2: vcf_2hap_merger_v2.py ---
# ############################################################################

def parse_info_step2(info):
    """
    Parse a semicolon-separated INFO string into a dictionary. (Step 2)
    """
    d = {}
    for item in info.strip().split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            d[key] = value
        else:
            d[item] = True
    return d

def format_info_step2(info_dict, order=None):
    """
    Convert an info dictionary back to a semicolon-separated string. (Step 2)
    """
    if order:
        items = []
        for key in order:
            if key in info_dict:
                items.append(f"{key}={info_dict[key]}")
        # Append any remaining keys
        for key, value in info_dict.items():
            if order and key not in order:
                items.append(f"{key}={value}")
        return ";".join(items)
    else:
        return ";".join(f"{k}={v}" for k, v in info_dict.items())

def read_vcf_step2(infile):
    """
    Read VCF file and return header lines and a list of entries. (Step 2)
    """
    header = []
    entries = []
    with open(infile, 'r') as f:
        for line in f:
            line = line.rstrip()
            if line.startswith("##"):
                header.append(line)
            elif line.startswith("#CHROM"):
                header.append(line)
            else:
                fields = line.split('\t')
                entry = {
                    "CHROM": fields[0],
                    "POS": int(fields[1]),
                    "ID": fields[2],
                    "REF": fields[3],
                    "ALT": fields[4],
                    "QUAL": fields[5],
                    "FILTER": fields[6],
                    "INFO": fields[7]
                }
                if len(fields) > 8:
                    entry["FORMAT"] = fields[8]
                    entry["SAMPLE"] = fields[9]
                else:
                    entry["FORMAT"] = "GT"
                    entry["SAMPLE"] = "."
                if entry["FILTER"] == "PASS_PRECISE":
                    entry["FILTER"] = "PRECISE"
                elif entry["FILTER"] == "PASS_IMPRECISE":
                    entry["FILTER"] = "IMPRECISE"
                entries.append(entry)
    return header, entries

def rename_hap_tags_step2(info, hap):
    """
    Rename haplotype-specific tags to include suffix _PAT or _MAT. (Step 2)
    """
    new_info = {}
    for key, value in info.items():
        if key in ["SUPSTRAND_END", "VALID_UNIQ_CLUSTER_END", "SUPMQ_END",
                   "SUPSTRAND_START", "VALID_UNIQ_CLUSTER_START", "SUPMQ_START"]:
            new_info[f"{key}_{hap}"] = value
        elif key in ["CIPOS", "CIEND"]:
            new_info[key] = value
        else:
            new_info[key] = value
    return new_info

def merge_group_step2(entries):
    """
    Given a list of entries for the same ORIGINAL_LOCUS, merge them. (Step 2)
    """
    if len(entries) == 1:
        entry = entries[0]
        info = parse_info_step2(entry["INFO"])
        source = info.get("SOURCE", "").upper()
        if "NOSIGNAL" in entry["FILTER"]:
            entry["SAMPLE"] = "0/0"
        if source == "PATERNAL":
            info = rename_hap_tags_step2(info, "PAT")
            info["CIPOS_PAT"] = info.get("CIPOS", "")
            info["CIEND_PAT"] = info.get("CIEND", "")
        elif source == "MATERNAL":
            info = rename_hap_tags_step2(info, "MAT")
            info["CIPOS_MAT"] = info.get("CIPOS", "")
            info["CIEND_MAT"] = info.get("CIEND", "")
        entry["INFO"] = format_info_step2(info)
        return [entry]
    
    elif len(entries) == 2:
        e1 = entries[0]
        e2 = entries[1]
        info1 = parse_info_step2(e1["INFO"])
        info2 = parse_info_step2(e2["INFO"])
        source1 = info1.get("SOURCE", "").upper()
        source2 = info2.get("SOURCE", "").upper()
        info1 = rename_hap_tags_step2(info1, "PAT" if source1 == "PATERNAL" else "MAT")
        info2 = rename_hap_tags_step2(info2, "PAT" if source2 == "PATERNAL" else "MAT")

        is_nosig_1 = ("NOSIGNAL" in e1["FILTER"])
        is_nosig_2 = ("NOSIGNAL" in e2["FILTER"])
        is_precise_1 = (e1["FILTER"] == "PRECISE")
        is_precise_2 = (e2["FILTER"] == "PRECISE")

        if is_nosig_1 and not is_nosig_2:
            if is_precise_2:
                e2_info = parse_info_step2(e2["INFO"])
                e2_info["MERGED_HAP"] = "PRECISE_CONF"
                e2["INFO"] = format_info_step2(e2_info)
                e2["SAMPLE"] = "0/1"
                return [e2]
            else:
                e2_info = parse_info_step2(e2["INFO"])
                e2_info["MERGED_HAP"] = "IMPRECISE_One_NoSignal"
                e2["INFO"] = format_info_step2(e2_info)
                e2["SAMPLE"] = "0/1"
                return [e2]
        elif is_nosig_2 and not is_nosig_1:
            if is_precise_1:
                e1_info = parse_info_step2(e1["INFO"])
                e1_info["MERGED_HAP"] = "PRECISE_CONF"
                e1["INFO"] = format_info_step2(e1_info)
                e1["SAMPLE"] = "0/1"
                return [e1]
            else:
                e1_info = parse_info_step2(e1["INFO"])
                e1_info["MERGED_HAP"] = "IMPRECISE_One_NoSignal"
                e1["INFO"] = format_info_step2(e1_info)
                e1["SAMPLE"] = "0/1"
                return [e1]
        elif is_nosig_1 and is_nosig_2:
            info_nosig = parse_info_step2(e1["INFO"])
            if "SOURCE" in info_nosig:
                del info_nosig["SOURCE"]
            e1["INFO"] = format_info_step2(info_nosig)
            e1["SAMPLE"] = "./."
            return [e1]
        
        cipos1 = info1.get("CIPOS", "")
        ciend1 = info1.get("CIEND", "")
        cipos2 = info2.get("CIPOS", "")
        ciend2 = info2.get("CIEND", "")
        start1 = e1["POS"]
        start2 = e2["POS"]
        try:
            end1 = int(info1.get("END", e1["POS"]))
        except:
            end1 = e1["POS"]
        try:
            end2 = int(info2.get("END", e2["POS"]))
        except:
            end2 = e2["POS"]

        if (("PRECISE" in e1["FILTER"] or "IMPRECISE" in e1["FILTER"]) and
            ("PRECISE" in e2["FILTER"] or "IMPRECISE" in e2["FILTER"])):
            
            if e1["FILTER"] == "IMPRECISE" or e2["FILTER"] == "IMPRECISE":
                call_type = "imprecise"
            else:
                call_type = "precise"
            diff_start = abs(start1 - start2)
            diff_end = abs(end1 - end2)
            
            if call_type == "precise":
                if diff_start == 0 and diff_end == 0:
                    merged_hap = "PRECISE_CONF"
                    merged_start = start1
                    merged_end = end1
                elif diff_start <= 500 and diff_end <= 500:
                    merged_hap = "PRECISE_Mid_CONF"
                    merged_start = int(round((start1 + start2) / 2))
                    merged_end = int(round((end1 + end2) / 2))
                else:
                    e1_info = parse_info_step2(e1["INFO"])
                    e2_info = parse_info_step2(e2["INFO"])
                    e1_info["MERGED_HAP"] = "PRECISE_Low_CONF"
                    e2_info["MERGED_HAP"] = "PRECISE_Low_CONF"
                    e1["INFO"] = format_info_step2(e1_info)
                    e2["INFO"] = format_info_step2(e2_info)
                    e1["SAMPLE"] = "0/1"
                    e2["SAMPLE"] = "0/1"
                    e1["FILTER"] = "PRECISE"
                    e2["FILTER"] = "PRECISE"
                    return [e1, e2]
            else:  # imprecise call
                if diff_start == 0 and diff_end == 0:
                    merged_hap = "IMPRECISE_Two_Hap"
                    merged_start = start1
                    merged_end = end1
                elif diff_start <= 500 and diff_end <= 500:
                    merged_hap = "IMPRECISE_Two_Hap"
                    merged_start = int(round((start1 + start2) / 2))
                    merged_end = int(round((end1 + end2) / 2))
                else:
                    e1["SAMPLE"] = "0/1"
                    e2["SAMPLE"] = "0/1"
                    e1["INFO"] = format_info_step2(info1)
                    e2["INFO"] = format_info_step2(info2)
                    e1["FILTER"] = "IMPRECISE"
                    e2["FILTER"] = "IMPRECISE"
                    return [e1, e2]

            merged_svlen = abs(merged_end - merged_start)
            merged_cipos_lower = min(start1, start2) - merged_start
            merged_cipos_upper = max(start1, start2) - merged_start
            merged_cipos = f"{merged_cipos_lower},{merged_cipos_upper}"
            merged_ciend_lower = min(end1, end2) - merged_end
            merged_ciend_upper = max(end1, end2) - merged_end
            merged_ciend = f"{merged_ciend_lower},{merged_ciend_upper}"

            info_merged = {}
            info_merged["END"] = str(merged_end)
            info_merged["SVLEN"] = str(merged_svlen)
            info_merged["SVTYPE"] = "INV"
            info_merged["MERGED_HAP"] = merged_hap
            info_merged["ORIGINAL_LOCUS"] = info1.get("ORIGINAL_LOCUS", info2.get("ORIGINAL_LOCUS", ""))
            info_merged["LABEL"] = info1.get("LABEL", info2.get("LABEL", "INV"))
            info_merged["CIPOS"] = merged_cipos
            info_merged["CIEND"] = merged_ciend
            info_merged["CIPOS_PAT"] = cipos1 if source1 == "PATERNAL" else cipos2
            info_merged["CIEND_PAT"] = ciend1 if source1 == "PATERNAL" else ciend2
            info_merged["SUPSTRAND_END_PAT"] = info1.get("SUPSTRAND_END_PAT", ".") if source1 == "PATERNAL" else info2.get("SUPSTRAND_END_PAT", ".")
            info_merged["VALID_UNIQ_CLUSTER_END_PAT"] = info1.get("VALID_UNIQ_CLUSTER_END_PAT", ".") if source1 == "PATERNAL" else info2.get("VALID_UNIQ_CLUSTER_END_PAT", ".")
            info_merged["SUPMQ_END_PAT"] = info1.get("SUPMQ_END_PAT", ".") if source1 == "PATERNAL" else info2.get("SUPMQ_END_PAT", ".")
            info_merged["SUPSTRAND_START_PAT"] = info1.get("SUPSTRAND_START_PAT", ".") if source1 == "PATERNAL" else info2.get("SUPSTRAND_START_PAT", ".")
            info_merged["VALID_UNIQ_CLUSTER_START_PAT"] = info1.get("VALID_UNIQ_CLUSTER_START_PAT", ".") if source1 == "PATERNAL" else info2.get("VALID_UNIQ_CLUSTER_START_PAT", ".")
            info_merged["SUPMQ_START_PAT"] = info1.get("SUPMQ_START_PAT", ".") if source1 == "PATERNAL" else info2.get("SUPMQ_START_PAT", ".")
            info_merged["CIPOS_MAT"] = cipos1 if source1 == "MATERNAL" else cipos2
            info_merged["CIEND_MAT"] = ciend1 if source1 == "MATERNAL" else ciend2
            info_merged["SUPSTRAND_END_MAT"] = info1.get("SUPSTRAND_END_MAT", ".") if source1 == "MATERNAL" else info2.get("SUPSTRAND_END_MAT", ".")
            info_merged["VALID_UNIQ_CLUSTER_END_MAT"] = info1.get("VALID_UNIQ_CLUSTER_END_MAT", ".") if source1 == "MATERNAL" else info2.get("VALID_UNIQ_CLUSTER_END_MAT", ".")
            info_merged["SUPMQ_END_MAT"] = info1.get("SUPMQ_END_MAT", ".") if source1 == "MATERNAL" else info2.get("SUPMQ_END_MAT", ".")
            info_merged["SUPSTRAND_START_MAT"] = info1.get("SUPSTRAND_START_MAT", ".") if source1 == "MATERNAL" else info2.get("SUPSTRAND_START_MAT", ".")
            info_merged["VALID_UNIQ_CLUSTER_START_MAT"] = info1.get("VALID_UNIQ_CLUSTER_START_MAT", ".") if source1 == "MATERNAL" else info2.get("VALID_UNIQ_CLUSTER_START_MAT", ".")
            info_merged["SUPMQ_START_MAT"] = info1.get("SUPMQ_START_MAT", ".") if source1 == "MATERNAL" else info2.get("SUPMQ_START_MAT", ".")

            order = ["END", "SVLEN", "SVTYPE", "MERGED_HAP", "ORIGINAL_LOCUS", "LABEL", "CIPOS", "CIEND",
                     "CIPOS_PAT", "CIEND_PAT", "SUPSTRAND_END_PAT", "VALID_UNIQ_CLUSTER_END_PAT", "SUPMQ_END_PAT",
                     "SUPSTRAND_START_PAT", "VALID_UNIQ_CLUSTER_START_PAT", "SUPMQ_START_PAT",
                     "CIPOS_MAT", "CIEND_MAT", "SUPSTRAND_END_MAT", "VALID_UNIQ_CLUSTER_END_MAT", "SUPMQ_END_MAT",
                     "SUPSTRAND_START_MAT", "VALID_UNIQ_CLUSTER_START_MAT", "SUPMQ_START_MAT"]
            merged_info_str = format_info_step2(info_merged, order)
            merged_entry = {
                "CHROM": e1["CHROM"],
                "POS": merged_start,
                "ID": ".",
                "REF": "N",
                "ALT": "<INV>",
                "QUAL": ".",
                "FILTER": e1["FILTER"],
                "INFO": merged_info_str,
                "FORMAT": e1.get("FORMAT", "GT"),
                "SAMPLE": "1/1"
            }
            return [merged_entry]
        else:
            gt1 = "0/1" if "NO_SIGNAL" not in e1["FILTER"] else "0/0"
            gt2 = "0/1" if "NO_SIGNAL" not in e2["FILTER"] else "0/0"
            e1["SAMPLE"] = gt1
            e2["SAMPLE"] = gt2
            e1["INFO"] = format_info_step2(info1)
            e2["INFO"] = format_info_step2(info2)
            return [e1, e2]
    else:
        return entries

def merge_vcf_step2(input_vcf, output_vcf):
    """
    Core logic for Step 2.
    """
    header, entries = read_vcf_step2(input_vcf)
    groups = {}
    for entry in entries:
        info = parse_info_step2(entry["INFO"])
        locus = info.get("ORIGINAL_LOCUS", f"{entry['CHROM']}:{entry['POS']}")
        groups.setdefault(locus, []).append(entry)
    merged_entries = []
    for locus, group in groups.items():
        merged_entries.extend(merge_group_step2(group))
    
    with open(output_vcf, "w") as out:
        default_header_lines = [
            "##fileformat=VCFv4.2",
            "##source=vcf_merger_v2.py",
            "##contig=<ID=chr1,length=248956422>",
            "##contig=<ID=chr2,length=242193529>",
            "##contig=<ID=chr3,length=198295559>",
            "##contig=<ID=chr4,length=190214555>",
            "##contig=<ID=chr5,length=181538259>",
            "##contig=<ID=chr6,length=170805979>",
            "##contig=<ID=chr7,length=159345973>",
            "##contig=<ID=chr8,length=145138636>",
            "##contig=<ID=chr9,length=138394717>",
            "##contig=<ID=chr10,length=133797422>",
            "##contig=<ID=chr11,length=135086622>",
            "##contig=<ID=chr12,length=133275309>",
            "##contig=<ID=chr13,length=114364328>",
            "##contig=<ID=chr14,length=107043718>",
            "##contig=<ID=chr15,length=101991189>",
            "##contig=<ID=chr16,length=90338345>",
            "##contig=<ID=chr17,length=83257441>",
            "##contig=<ID=chr18,length=80373285>",
            "##contig=<ID=chr19,length=58617616>",
            "##contig=<ID=chr20,length=64444167>",
            "##contig=<ID=chr21,length=46709983>",
            "##contig=<ID=chr22,length=50818468>",
            "##contig=<ID=chrX,length=156040895>",
            "##contig=<ID=chrY,length=57227415>",
            "##contig=<ID=chrM,length=16569>",
            "##FILTER=<ID=PRECISE,Description=\"Variant call is precise\">",
            "##FILTER=<ID=IMPRECISE,Description=\"Variant call is imprecise\">",
            "##FILTER=<ID=NOSIGNAL,Description=\"No haplotype signal detected\">",
            "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant\">",
            "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of the structural variant\">",
            "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">",
            "##INFO=<ID=ORIGINAL_LOCUS,Number=1,Type=String,Description=\"Original locus of the variant\">",
            "##INFO=<ID=LABEL,Number=1,Type=String,Description=\"Label for the variant\">",
            "##INFO=<ID=CIPOS,Number=2,Type=Integer,Description=\"Confidence interval around reported POS\">",
            "##INFO=<ID=CIEND,Number=2,Type=Integer,Description=\"Confidence interval around reported END\">",
            "##INFO=<ID=SUPMQ_START_MAT,Number=1,Type=String,Description=\"Support mapping quality at haplotype maternal start\">",
            "##INFO=<ID=SUPSTRAND_START_MAT,Number=1,Type=String,Description=\"Support strand at haplotype maternal start\">",
            "##INFO=<ID=SUPMQ_END_MAT,Number=1,Type=String,Description=\"Support mapping quality at haplotype maternal end\">",
            "##INFO=<ID=SUPSTRAND_END_MAT,Number=1,Type=String,Description=\"Support strand at haplotype maternal end\">",
            "##INFO=<ID=VALID_UNIQ_CLUSTER_START_MAT,Number=1,Type=String,Description=\"Valid unique cluster at haplotype maternal start\">",
            "##INFO=<ID=VALID_UNIQ_CLUSTER_END_MAT,Number=1,Type=String,Description=\"Valid unique cluster at haplotype maternal end\">",
            "##INFO=<ID=SUPMQ_START_PAT,Number=1,Type=String,Description=\"Support mapping quality at haplotype paternal start\">",
            "##INFO=<ID=SUPSTRAND_START_PAT,Number=1,Type=String,Description=\"Support strand at haplotype paternal start\">",
            "##INFO=<ID=SUPMQ_END_PAT,Number=1,Type=String,Description=\"Support mapping quality at haplotype paternal end\">",
            "##INFO=<ID=SUPSTRAND_END_PAT,Number=1,Type=String,Description=\"Support strand at haplotype paternal end\">",
            "##INFO=<ID=VALID_UNIQ_CLUSTER_START_PAT,Number=1,Type=String,Description=\"Valid unique cluster at haplotype paternal start\">",
            "##INFO=<ID=VALID_UNIQ_CLUSTER_END_PAT,Number=1,Type=String,Description=\"Valid unique cluster at haplotype paternal end\">",
            "##INFO=<ID=CIPOS_PAT,Number=1,Type=String,Description=\"Confidence interval around start position, paternal\">",
            "##INFO=<ID=CIEND_PAT,Number=1,Type=String,Description=\"Confidence interval around end position, paternal\">",
            "##INFO=<ID=CIPOS_MAT,Number=1,Type=String,Description=\"Confidence interval around start position, maternal\">",
            "##INFO=<ID=CIEND_MAT,Number=1,Type=String,Description=\"Confidence interval around end position, maternal\">",
            "##INFO=<ID=VALID_CLUSTER_PAIR,Number=1,Type=String,Description=\"Valid cluster pair information\">",
            "##INFO=<ID=MERGED_HAP,Number=1,Type=String,Description=\"Merged haplotype call type\">",
            "##INFO=<ID=SOURCE,Number=1,Type=String,Description=\"Source haplotype information. Used for heterozygous variants where one haplotype has signal and the other is NOSIGNAL\">",
            "##INFO=<ID=SUPMQ_START,Number=1,Type=String,Description=\"Support mapping quality at haplotype start for heterozygous variants (signal from one haplotype)\">",
            "##INFO=<ID=SUPSTRAND_START,Number=1,Type=String,Description=\"Support strand at haplotype start for heterozygous variants (signal from one haplotype)\">",
            "##INFO=<ID=SUPMQ_END,Number=1,Type=String,Description=\"Support mapping quality at haplotype end for heterozygous variants (signal from one haplotype)\">",
            "##INFO=<ID=SUPSTRAND_END,Number=1,Type=String,Description=\"Support strand at haplotype end for heterozygous variants (signal from one haplotype)\">",
            "##INFO=<ID=VALID_UNIQ_CLUSTER_START,Number=1,Type=String,Description=\"Valid unique cluster at haplotype start for heterozygous variants (signal from one haplotype)\">",
            "##INFO=<ID=VALID_UNIQ_CLUSTER_END,Number=1,Type=String,Description=\"Valid unique cluster at haplotype end for heterozygous variants (signal from one haplotype)\">",
            "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"
        ]
        
        for line in default_header_lines:
            out.write(line + "\n")

        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")

        for entry in merged_entries:
            fields = [
                entry["CHROM"],
                str(entry["POS"]),
                entry["ID"],
                entry["REF"],
                entry["ALT"],
                entry["QUAL"],
                entry["FILTER"],
                entry["INFO"],
                entry.get("FORMAT", "GT"),
                entry.get("SAMPLE", ".")
            ]
            out.write("\t".join(fields) + "\n")

def run_step2(input_vcf, output_vcf):
    """
    Runner function for Step 2.
    """
    try:
        merge_vcf_step2(input_vcf, output_vcf)
    except Exception as e:
        print(f"Error in Step 2 (merge_vcf) for {input_vcf}: {e}", file=sys.stderr)
        sys.exit(1)


# ############################################################################
# --- STEP 3: rm_redundent_inv.py ---
# ############################################################################

# Define a named tuple for BedEntry for Step 3
BedEntry_step3 = collections.namedtuple('BedEntry_step3', ['chrom', 'start', 'end', 'original_line'])

def parse_bed_line_step3(line_number, line):
    """
    Parses a line from a BED file (Step 3).
    """
    parts = line.strip().split()
    if len(parts) < 3:
        print(f"Warning (Step 3): Skipping malformed BED line {line_number}: '{line.strip()}'", file=sys.stderr)
        return None
    try:
        chrom = parts[0]
        start = int(parts[1])
        end = int(parts[2])
        return BedEntry_step3(chrom, start, end, line.strip())
    except ValueError:
        print(f"Warning (Step 3): Skipping malformed BED line {line_number}: '{line.strip()}'", file=sys.stderr)
        return None

def generate_filter_pattern_from_bed_step3(inversions_bed_filepath, log_filepath):
    """
    Processes the BED file to identify entries strictly contained within others (Step 3).
    """
    all_entries_by_chrom = collections.defaultdict(list)
    ordered_chroms_from_input = [] 
    chry_filtered_log_entries = []
    
    bed_processing_log_header = f"--- (Step 3) Processing BED file: {inversions_bed_filepath} ---\n"
    
    try:
        with open(inversions_bed_filepath, 'r') as f_in:
            for i, line in enumerate(f_in):
                line_num = i + 1
                entry = parse_bed_line_step3(line_num, line)
                if not entry:
                    continue

                if entry.chrom == 'chrY':
                    chry_filtered_log_entries.append(f"Filtered out chrY entry from BED: {entry.original_line}")
                    continue

                if entry.chrom not in all_entries_by_chrom:
                    if entry.chrom not in ordered_chroms_from_input:
                         ordered_chroms_from_input.append(entry.chrom)
                
                all_entries_by_chrom[entry.chrom].append(entry)

    except FileNotFoundError:
        error_msg = f"Error (Step 3): Input BED file not found at {inversions_bed_filepath}"
        print(error_msg, file=sys.stderr)
        try:
            with open(log_filepath, 'w') as f_log:
                f_log.write(bed_processing_log_header)
                f_log.write(error_msg + "\n")
        except IOError:
            print(f"Critical Error (Step 3): Could not write initial error to log file {log_filepath}", file=sys.stderr)
        return None

    selected_entry_identifiers = set() 

    for chrom in ordered_chroms_from_input: 
        entries_for_chrom = all_entries_by_chrom.get(chrom, [])
        if not entries_for_chrom or len(entries_for_chrom) < 2:
            continue

        for i in range(len(entries_for_chrom)):
            e1 = entries_for_chrom[i]
            e1_identifier = f"{e1.chrom}_{e1.start}_{e1.end}"

            for j in range(len(entries_for_chrom)):
                if i == j: 
                    continue
                
                e2 = entries_for_chrom[j]
                overlap = (max(e1.start, e2.start) < min(e1.end, e2.end))

                if overlap:
                    e1_is_contained_in_e2 = (e2.start <= e1.start and e1.end <= e2.end)
                    if e1_is_contained_in_e2:
                        e1_is_strictly_contained_in_e2 = (e1.start > e2.start or e1.end < e2.end)
                        if e1_is_strictly_contained_in_e2:
                            selected_entry_identifiers.add(e1_identifier)
    
    final_result_parts = []
    parsed_selected_entries = collections.defaultdict(list)
    for entry_str in selected_entry_identifiers:
        parts = entry_str.split('_')
        if len(parts) >= 3:
            chrom_name = parts[0] 
            try:
                start_pos = int(parts[1])
                end_pos = int(parts[2])
                parsed_selected_entries[chrom_name].append((start_pos, end_pos, entry_str))
            except ValueError:
                print(f"Warning (Step 3): Could not parse identifier '{entry_str}'", file=sys.stderr)
        else:
            print(f"Warning (Step 3): Malformed identifier '{entry_str}'", file=sys.stderr)

    for chrom_name in ordered_chroms_from_input:
        if chrom_name in parsed_selected_entries:
            sorted_entries_for_chrom = sorted(parsed_selected_entries[chrom_name]) 
            for _, _, entry_str in sorted_entries_for_chrom:
                final_result_parts.append(entry_str)
                
    filter_pattern_string = "|".join(final_result_parts)

    try:
        with open(log_filepath, 'w') as f_log:
            f_log.write(bed_processing_log_header)
            if chry_filtered_log_entries:
                for log_entry in chry_filtered_log_entries:
                    f_log.write(log_entry + "\n")
            else:
                f_log.write("No chrY entries were filtered from BED.\n")
            f_log.write(f"Generated filter pattern from BED: {filter_pattern_string if filter_pattern_string else 'NONE'}\n")
            f_log.write("--- End of (Step 3) BED processing ---\n")
        print(f"(Step 3) BED processing log written to: {log_filepath}", file=sys.stderr)
    except IOError:
        print(f"Error (Step 3): Could not write BED processing details to log file {log_filepath}", file=sys.stderr)

    return filter_pattern_string

def filter_vcf_with_pattern_step3(input_vcf_filepath, output_vcf_filepath, filter_pattern_string, log_file_object):
    """
    Filters the input VCF file using the provided pattern (Step 3).
    """
    if not filter_pattern_string:
        message = f"No filter pattern provided for {input_vcf_filepath}. Copying VCF file as is to {output_vcf_filepath}."
        print(message, file=sys.stderr)
        log_file_object.write(message + "\n")
        try:
            with open(input_vcf_filepath, 'r') as f_in, open(output_vcf_filepath, 'w') as f_out:
                for line in f_in:
                    f_out.write(line)
            copy_summary = f"VCF copying for {input_vcf_filepath} complete (no pattern applied)."
            print(copy_summary, file=sys.stderr)
            log_file_object.write(copy_summary + "\n")
        except IOError as e:
            error_copy = f"Error copying VCF file {input_vcf_filepath}: {e}"
            print(error_copy, file=sys.stderr)
            log_file_object.write(error_copy + "\n")
        return

    try:
        regex_pattern = re.compile(filter_pattern_string)
    except re.error as e:
        regex_error_msg = f"Error compiling regex pattern for {input_vcf_filepath}: {e}. Copying file."
        print(regex_error_msg, file=sys.stderr)
        log_file_object.write(regex_error_msg + "\n")
        try:
            with open(input_vcf_filepath, 'r') as f_in, open(output_vcf_filepath, 'w') as f_out:
                for line in f_in:
                    f_out.write(line)
        except IOError as e_copy:
            error_fallback = f"Error copying VCF file {input_vcf_filepath} during fallback: {e_copy}"
            print(error_fallback, file=sys.stderr)
            log_file_object.write(error_fallback + "\n")
        return

    lines_written = 0
    lines_filtered = 0

    try:
        with open(input_vcf_filepath, 'r') as f_in, open(output_vcf_filepath, 'w') as f_out:
            for line in f_in:
                if line.startswith('#'): 
                    f_out.write(line)
                    lines_written +=1
                    continue
                
                if regex_pattern.search(line):
                    lines_filtered += 1
                    continue 
                
                f_out.write(line)
                lines_written +=1
        
        filter_summary = f"VCF filtering for {input_vcf_filepath} complete. Output: {output_vcf_filepath}. Lines written: {lines_written}, Lines filtered out: {lines_filtered}"
        print(filter_summary, file=sys.stderr)
        log_file_object.write(filter_summary + "\n")

    except FileNotFoundError:
        fnf_error = f"Error: Input VCF file not found at {input_vcf_filepath}"
        print(fnf_error, file=sys.stderr)
        log_file_object.write(fnf_error + "\n")
    except IOError as e:
        io_error = f"Error during VCF filtering for {input_vcf_filepath}: {e}"
        print(io_error, file=sys.stderr)
        log_file_object.write(io_error + "\n")

def run_step3(inversions_bed, input_vcf1, input_vcf2, log_file, out_vcf1_path, out_vcf2_path):
    """
    Runner function for Step 3.
    """
    print(f"(Step 3) Generating filter pattern from BED file: {inversions_bed}", file=sys.stderr)
    filter_pattern = generate_filter_pattern_from_bed_step3(inversions_bed, log_file)

    if filter_pattern is None:
        print("Critical error during Step 3 BED processing. Aborting VCF filtering.", file=sys.stderr)
        sys.exit(1)
    
    print(f"(Step 3) Filter pattern generated: '{filter_pattern if filter_pattern else 'NONE'}'", file=sys.stderr)

    try:
        with open(log_file, 'a') as f_log:
            vcf_files_to_process = [
                (input_vcf1, out_vcf1_path), 
                (input_vcf2, out_vcf2_path)
            ]
            for vcf_input_path, vcf_output_path in vcf_files_to_process:
                log_header_vcf = f"\n--- (Step 3) Starting VCF processing for: {vcf_input_path} ---"
                print(log_header_vcf, file=sys.stderr)
                f_log.write(log_header_vcf + "\n")
                
                filter_vcf_with_pattern_step3(vcf_input_path, vcf_output_path, filter_pattern, f_log)
                
                log_footer_vcf = f"--- (Step 3) Finished VCF processing for: {vcf_input_path}. Output: {vcf_output_path} ---"
                print(log_footer_vcf, file=sys.stderr)
                f_log.write(log_footer_vcf + "\n\n")
                
    except IOError:
        print(f"Error (Step 3): Could not open or append to log file {log_file} for VCF processing.", file=sys.stderr)
        sys.exit(1)
        
    print("(Step 3) All processing complete.", file=sys.stderr)


# ############################################################################
# --- STEP 4: rp_centromere_location_annotate.py ---
# ############################################################################

def load_bed_step4_intervaltree(bed_filename):
    """
    Loads a BED file into an interval tree. (Step 4)
    """
    trees = {}
    try:
        with open(bed_filename, "r") as bed_file:
            for line in bed_file:
                if line.startswith("#") or line.strip() == "":
                    continue
                fields = line.strip().split()
                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])
                if chrom not in trees:
                    trees[chrom] = IntervalTree()
                trees[chrom].addi(start, end)
    except FileNotFoundError:
        print(f"Error (Step 4): Annotation BED file not found: {bed_filename}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error (Step 4): Could not parse BED file {bed_filename}: {e}", file=sys.stderr)
        sys.exit(1)
    return trees

def annotate_variant_step4(chrom, start0, end1, tree_dict):
    """
    Checks for overlap with an interval tree. (Step 4)
    """
    start_overlap = False
    end_overlap = False

    if chrom in tree_dict:
        tree = tree_dict[chrom]
        if tree.overlap(start0, start0 + 1):
            start_overlap = True
        if tree.overlap(end1 - 1, end1):
            end_overlap = True

    if start_overlap and not end_overlap:
        return "rp_left"
    elif end_overlap and not start_overlap:
        return "rp_right"
    elif start_overlap and end_overlap:
        return "rp2sides"
    else:
        return "."

def ensure_samples_in_header_step4(vcf_in_path):
    """
    Checks if VCF has sample columns. (Step 4)
    """
    try:
        with pysam.VariantFile(vcf_in_path) as vcf:
            if list(vcf.header.samples):
                return vcf
            else:
                # This script relies on the VCFs from Step 2, which SHOULD have samples.
                # If not, we can't use pysam properly.
                print(f"Warning (Step 4): VCF {vcf_in_path} has no sample columns. Re-opening without sample checks.", file=sys.stderr)
                # Fallback to text parsing? No, the original script uses pysam.
                # Let's check if the *original* script fails.
                # The original script *fails* if no samples. This is intended behavior.
                raise ValueError(f"VCF {vcf_in_path} does not contain sample genotype columns.")
    except Exception as e:
        print(f"Error (Step 4): Could not open VCF {vcf_in_path} with pysam. {e}", file=sys.stderr)
        sys.exit(1)

def run_step4(input_vcf, output_vcf, bed_5kb_path, bed_10kb_path, bed_centromeres_path):
    """
    Runner function for Step 4.
    """
    try:
        # Load BED files
        print(f"(Step 4) Loading annotation BEDs...")
        tree_5kb = load_bed_step4_intervaltree(bed_5kb_path)
        tree_10kb = load_bed_step4_intervaltree(bed_10kb_path)
        tree_centromere = load_bed_step4_intervaltree(bed_centromeres_path)
        print(f"(Step 4) Annotation BEDs loaded.")

        # Open input and output VCFs
        # Note: The original ensure_samples_in_header fails if no samples.
        # The VCF from step 2 *should* have a sample.
        vcf_in = pysam.VariantFile(input_vcf)
        
        header = vcf_in.header.copy()
        header.info.add("Repeat_ge5Kb", 1, "String", "Repeat overlap (>=5kb): rp_left, rp_right, rp2sides")
        header.info.add("Repeat_ge10Kb", 1, "String", "Repeat overlap (>=10kb): rp_left, rp_right, rp2sides")
        header.info.add("Centromere", 0, "Flag", "Variant overlaps centromere region")

        vcf_out = pysam.VariantFile(output_vcf, "w", header=header)

        for rec in vcf_in:
            chrom = rec.chrom
            start0 = rec.start
            try:
                end1 = rec.stop
            except KeyError:
                # Fallback if END tag is missing (though it shouldn't be from step 2)
                print(f"Warning (Step 4): Record {rec.chrom}:{rec.pos} missing END tag. Using POS+1.", file=sys.stderr)
                end1 = rec.pos + 1


            ann_5kb = annotate_variant_step4(chrom, start0, end1, tree_5kb)
            ann_10kb = annotate_variant_step4(chrom, start0, end1, tree_10kb)
            centromere_hit = annotate_variant_step4(chrom, start0, end1, tree_centromere) != "."

            if ann_5kb != ".":
                rec.info["Repeat_ge5Kb"] = ann_5kb
            if ann_10kb != ".":
                rec.info["Repeat_ge10Kb"] = ann_10kb
            if centromere_hit:
                rec.info["Centromere"] = True

            vcf_out.write(rec)

        vcf_out.close()
        vcf_in.close()
    
    except Exception as e:
        print(f"Error in Step 4 (annotate_vcf) for {input_vcf}: {e}", file=sys.stderr)
        sys.exit(1)


# ############################################################################
# --- STEP 5: merge2aligner.py ---
# ############################################################################

FAIL_HEADER_LINE_STEP5 = '##FILTER=<ID=FAIL,Description="Both callers NOSIGNAL; flagged as FAIL">'
KEEP_FILTER_STEP5      = "NOSIGNAL"
DROP_CHR_STEP5         = "chrY"
CENTROMERE_TAG_STEP5   = "Centromere"


def ensure_fail_filter_step5(hdr: pysam.VariantHeader):
    if "FAIL" not in hdr.filters:
        hdr.add_line(FAIL_HEADER_LINE_STEP5)

def first_filter_step5(rec: pysam.VariantRecord) -> str:
    return next(iter(rec.filter), "PASS")

def sub_category_step5(rec: pysam.VariantRecord) -> str:
    return rec.info.get("MERGED_HAP", first_filter_step5(rec))

def chrom_sort_key_step5(chrom: str) -> int:
    """Numeric sort for chr1..chr22, chrX, chrY, others last."""
    c = chrom.lower().removeprefix("chr")
    if c.isdigit():
        return int(c)
    if c == "x":
        return 23
    if c == "y":
        return 24
    return 25  # others

def record_sort_key_step5(rec: pysam.VariantRecord):
    return (chrom_sort_key_step5(rec.chrom), rec.pos)

def write_log_step5(path, mm2_kept, vac_kept, dropped_centromere, dropped_chrY, failed_loci):
    with open(path, "w") as lg:
        lg.write("merge2aligners2finalVCF log (Step 5)\n\n")

        lg.write(f"Dropped due to {CENTROMERE_TAG_STEP5}: {len(dropped_centromere)} loci\n")
        if dropped_centromere:
            lg.write("    loci: " + ", ".join(sorted(dropped_centromere)) + "\n")
        lg.write("\n")

        lg.write(f"Dropped on {DROP_CHR_STEP5}: {len(dropped_chrY)} loci\n")
        if dropped_chrY:
            lg.write("    loci: " + ", ".join(sorted(dropped_chrY)) + "\n")
        lg.write("\n")

        lg.write(f"Both NOSIGNAL → dropped: {len(failed_loci)} loci\n")
        if failed_loci:
            lg.write("    loci: " + ", ".join(sorted(failed_loci)) + "\n")
        lg.write("\n")

        total_entries_mm2 = sum(len(v) for v in mm2_kept.values())
        lg.write(f"Kept from minimap2: {total_entries_mm2} entries\n")
        for subcat, loci in sorted(mm2_kept.items()):
            distinct = sorted(set(loci))
            lg.write(f"{subcat}: {len(distinct)} loci\n")
            lg.write("    loci: " + ", ".join(distinct) + "\n")
            lg.write(f"    entries: {len(loci)}\n")
        lg.write("\n")

        total_entries_vac = sum(len(v) for v in vac_kept.values())
        lg.write(f"Kept from VACmap: {total_entries_vac} entries\n")
        for subcat, loci in sorted(vac_kept.items()):
            distinct = sorted(set(loci))
            lg.write(f"{subcat}: {len(distinct)} loci\n")
            lg.write("    loci: " + ", ".join(distinct) + "\n")
            lg.write(f"    entries: {len(loci)}\n")

def run_step5(minimap2_vcf, vacmap_vcf, merged_vcf):
    """
    Runner function for Step 5.
    """
    try:
        mm2 = pysam.VariantFile(minimap2_vcf)
        vac = pysam.VariantFile(vacmap_vcf)
        for hdr in (mm2.header, vac.header):
            ensure_fail_filter_step5(hdr)

        header_main    = mm2.header.copy()
        header_dropped = mm2.header.copy()
        main_recs      = []
        dropped_recs   = []

        mm2_kept            = defaultdict(list)
        vac_kept            = defaultdict(list)
        dropped_centromere  = []
        dropped_chrY        = []
        failed_loci         = []

        mm2_by_locus = defaultdict(list)
        for rec in mm2:
            mm2_by_locus[rec.info["ORIGINAL_LOCUS"]].append(rec)

        vac_by_locus = defaultdict(list)
        for rec in vac:
            vac_by_locus[rec.info["ORIGINAL_LOCUS"]].append(rec)

        all_loci = sorted(set(mm2_by_locus) | set(vac_by_locus))
        for locus in all_loci:
            recs_mm2 = mm2_by_locus.get(locus, [])
            recs_vac = vac_by_locus.get(locus, [])

            combined = recs_mm2 + recs_vac
            if not combined:
                continue

            # 1) Centromere filter
            if any(CENTROMERE_TAG_STEP5 in rec.info for rec in combined):
                dropped_centromere.append(locus)
                dropped_recs.extend(combined)
                continue

            # 2) chrY filter
            rec0 = recs_mm2[0] if recs_mm2 else recs_vac[0]
            if rec0.chrom == DROP_CHR_STEP5:
                dropped_chrY.append(locus)
                dropped_recs.extend(combined)
                continue

            # 3) If any minimap2 record != NOSIGNAL -> keep all minimap2
            mm2_filt = [first_filter_step5(r) for r in recs_mm2]
            if recs_mm2 and any(f != KEEP_FILTER_STEP5 for f in mm2_filt):
                for r, f in zip(recs_mm2, mm2_filt):
                    main_recs.append(r)
                    mm2_kept[sub_category_step5(r)].append(locus)
                continue

            # 4) Otherwise fallback to VACmap
            vac_filt = [first_filter_step5(r) for r in recs_vac]
            if recs_vac and any(f != KEEP_FILTER_STEP5 for f in vac_filt):
                for r, f in zip(recs_vac, vac_filt):
                    main_recs.append(r)
                    vac_kept[sub_category_step5(r)].append(locus)
            elif recs_vac:
                # both NOSIGNAL -> drop them
                failed_loci.append(locus)
                dropped_recs.extend(combined)

        main_recs_sorted    = sorted(main_recs,    key=record_sort_key_step5)
        dropped_recs_sorted = sorted(dropped_recs, key=record_sort_key_step5)

        out_main = pysam.VariantFile(merged_vcf, "w", header=header_main)
        for rec in main_recs_sorted:
            out_main.write(rec)
        out_main.close()

        dropped_vcf = merged_vcf + ".dropped.vcf"
        out_drop = pysam.VariantFile(dropped_vcf, "w", header=header_dropped)
        for rec in dropped_recs_sorted:
            out_drop.write(rec)
        out_drop.close()
        
        log_path = merged_vcf + ".log"
        write_log_step5(
            log_path,
            mm2_kept, vac_kept,
            dropped_centromere, dropped_chrY, failed_loci
        )

        print(f"✔ (Step 5) main output:    {merged_vcf}")
        print(f"✔ (Step 5) dropped output: {dropped_vcf}")
        print(f"✔ (Step 5) log:            {log_path}")
    
    except Exception as e:
        print(f"Error in Step 5 (merge2aligner) for {minimap2_vcf} / {vacmap_vcf}: {e}", file=sys.stderr)
        sys.exit(1)


# ############################################################################
# --- STEP 6: stat_mergedVCF_v2.py ---
# ############################################################################

def parse_info_step6(info):
    """
    Parse a semicolon-separated INFO string into a dictionary. (Step 6)
    """
    d = {}
    for item in info.strip().split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            d[key] = value
        else:
            d[item] = True
    return d

def safe_parse_pair_step6(value):
    """
    Given "a,b", return (float(a), float(b)). (Step 6)
    """
    try:
        parts = value.split(',')
        return (float(parts[0]), float(parts[1]))
    except Exception:
        return None

def average_pair_step6(pairs):
    """
    Given [(a, b), ...], return (avg_a, avg_b). (Step 6)
    """
    if not pairs:
        return (None, None)
    sum_a = sum(a for a, b in pairs)
    sum_b = sum(b for a, b in pairs)
    n = len(pairs)
    return (sum_a/n, sum_b/n)

def median_pair_step6(pairs):
    """
    Given [(a, b), ...], return (med_a, med_b). (Step 6)
    """
    if not pairs:
        return (None, None)
    lowers = [a for a, b in pairs]
    uppers = [b for a, b in pairs]
    return (statistics.median(lowers), statistics.median(uppers))

def process_vcf_step6(infile):
    """
    Parse the VCF file and collect statistics. (Step 6)
    """
    stats = {}
    repeats = {}
    total = 0

    length_bins_def = [
        (50, 500, "50-500 bp"),
        (500, 1000, "500 bp - 1 Kb"),
        (1000, 10000, "1-10 Kb"),
        (10000, 50000, "10-50 Kb"),
        (50000, 100000, "50-100 Kb"),
        (100000, 200000, "100-200 Kb"),
        (200000, 500000, "200-500 Kb"),
        (500000, 1000000, "500 Kb - 1 Mb"),
        (1000000, 2000000, "1-2 Mb"),
        (2000000, 5000000, "2-5 Mb"),
        (5000000, float("inf"), ">5 Mb")
    ]

    try:
        with open(infile, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 8:
                    fields = line.split()
                if len(fields) < 8:
                    continue

                if len(fields) >= 10:
                    chrom, pos, vid, ref, alt, qual, fltr, info_str, fmt, sample = fields[:10]
                else:
                    chrom, pos, vid, ref, alt, qual, fltr, info_str = fields[:8]
                    sample = "NA"

                try:
                    pos = int(pos)
                except Exception:
                    continue

                info = parse_info_step6(info_str)
                total += 1

                if sample == "./.":
                    category = "NO_SIGNAL"
                elif "MERGED_HAP" in info:
                    category = info["MERGED_HAP"]
                else:
                    category = fltr

                if category not in stats:
                    stats[category] = {
                        "count": 0,
                        "cipos": [],
                        "ciend": [],
                        "gt_counts": {},
                        "svlen": [],
                        "length_bins": { bin_label: 0 for (_, _, bin_label) in length_bins_def },
                        "repeat_ge10Kb": 0,
                        "repeat_ge5Kb": 0,
                        "centromere": 0
                    }
                stats[category]["count"] += 1

                orig = info.get("ORIGINAL_LOCUS", "NA")
                if category not in repeats:
                    repeats[category] = {}
                repeats[category][orig] = repeats[category].get(orig, 0) + 1

                cipos_val = info.get("CIPOS", None)
                if cipos_val:
                    pair = safe_parse_pair_step6(cipos_val)
                    if pair:
                        stats[category]["cipos"].append(pair)
                
                ciend_val = info.get("CIEND", None)
                if ciend_val:
                    pair = safe_parse_pair_step6(ciend_val)
                    if pair:
                        stats[category]["ciend"].append(pair)
                
                svlen_val = info.get("SVLEN", None)
                if svlen_val:
                    try:
                        svlen_num = abs(float(svlen_val))
                        stats[category]["svlen"].append(svlen_num)
                        for lower, upper, label in length_bins_def:
                            if lower <= svlen_num < upper:
                                stats[category]["length_bins"][label] += 1
                                break
                    except Exception:
                        pass
                
                if 'Repeat_ge10Kb' in info:
                    stats[category]["repeat_ge10Kb"] += 1
                if 'Repeat_ge5Kb' in info:
                    stats[category]["repeat_ge5Kb"] += 1
                if 'Centromere' in info:
                    stats[category]["centromere"] += 1
                
                stats[category]["gt_counts"][sample] = stats[category]["gt_counts"].get(sample, 0) + 1
    
    except FileNotFoundError:
        print(f"Error (Step 6): Input VCF not found: {infile}", file=sys.stderr)
        return {}, 0, {}
    except Exception as e:
        print(f"Error (Step 6) processing VCF {infile}: {e}", file=sys.stderr)
        return {}, 0, {}

    return stats, total, repeats

def print_stats_step6(stats, total, repeats):
    """
    Prints the stats. (Step 6)
    """
    print("VCF Statistics Summary (Step 6)")
    print("======================")
    print(f"Total number of INV records: {total}\n")
    for category, data in stats.items():
        count = data["count"]
        prop = (count / total) * 100 if total > 0 else 0
        avg_cipos = average_pair_step6(data["cipos"])
        med_cipos = median_pair_step6(data["cipos"])
        avg_ciend = average_pair_step6(data["ciend"])
        med_ciend = median_pair_step6(data["ciend"])
        print(f"Category: {category}")
        print(f"\tCount: {count} ({prop:.2f}% of total)")
        if category in repeats:
            loci_counts = repeats[category]
            unique_loci = len(loci_counts)
            total_occurrences = sum(loci_counts.values())
            duplicate_occurrences = total_occurrences - unique_loci
            print(f"\tRepeat info: Unique INV events: {unique_loci}; Duplicate occurrences: {duplicate_occurrences}")
        else:
            print("\tRepeat info: N/A")
        if avg_cipos[0] is not None:
            print(f"\tAverage CIPOS: lower={avg_cipos[0]:.2f}, upper={avg_cipos[1]:.2f}")
            print(f"\tMedian CIPOS: lower={med_cipos[0]:.2f}, upper={med_cipos[1]:.2f}")
        else:
            print("\tAverage CIPOS: N/A")
            print("\tMedian CIPOS: N/A")
        if avg_ciend[0] is not None:
            print(f"\tAverage CIEND: lower={avg_ciend[0]:.2f}, upper={avg_ciend[1]:.2f}")
            print(f"\tMedian CIEND: lower={med_ciend[0]:.2f}, upper={med_ciend[1]:.2f}")
        else:
            print("\tAverage CIEND: N/A")
            print("\tMedian CIEND: N/A")
        print("\tGenotype counts:")
        total_gt = sum(data["gt_counts"].values())
        for gt, cnt in data["gt_counts"].items():
            perc = (cnt / total_gt) * 100 if total_gt > 0 else 0
            print(f"\t\t{gt}: {cnt} ({perc:.2f}%)")
        if data["svlen"]:
            avg_svlen = sum(data["svlen"]) / len(data["svlen"])
            print(f"\tAverage SVLEN: {avg_svlen:.2f}")
            print("\tSVLEN distribution by length bins:")
            for bin_label, cnt in data["length_bins"].items():
                print(f"\t\t{bin_label}: {cnt}")
        else:
            print("\tAverage SVLEN: N/A")
            print("\tSVLEN distribution by length bins: N/A")
        
        print(f"\tRepeat_ge10Kb: {data['repeat_ge10Kb']}")
        print(f"\tRepeat_ge5Kb: {data['repeat_ge5Kb']}")
        print(f"\tCentromere: {data['centromere']}")
        print("")

def run_step6(input_vcf, output_stats_file):
    """
    Runner function for Step 6.
    Redirects stdout to a file.
    """
    stats, total, repeats = process_vcf_step6(input_vcf)
    
    if total == 0:
        print(f"(Step 6) No records found in {input_vcf}. Skipping stats generation.", file=sys.stderr)
        return

    original_stdout = sys.stdout
    try:
        with open(output_stats_file, 'w') as f:
            sys.stdout = f
            print_stats_step6(stats, total, repeats)
    except Exception as e:
        print(f"Error (Step 6) writing stats to {output_stats_file}: {e}", file=sys.stderr)
    finally:
        sys.stdout = original_stdout # Restore stdout


# ############################################################################
# --- STEP 7: table_matrix_v2.py ---
# ############################################################################

def parse_info_field_step7(info_str):
    """
    Parses the INFO field of a VCF record into a dictionary. (Step 7)
    """
    info_dict = {}
    if info_str == ".":
        return info_dict
    fields = info_str.split(';')
    for field in fields:
        if '=' in field:
            key, value = field.split('=', 1)
            info_dict[key] = value
        else:
            info_dict[field] = True
    return info_dict

def determine_category_step7(filter_val, info_dict):
    """
    Determines the category of a VCF record. (Step 7)
    """
    if filter_val == "NOSIGNAL":
        return "NOSIGNAL"
    elif filter_val == "IMPRECISE":
        merged_hap = info_dict.get("MERGED_HAP")
        if merged_hap == "IMPRECISE_Two_Hap":
            return "IMPRECISE_Two_Hap"
        elif merged_hap == "IMPRECISE_One_NoSignal":
            return "IMPRECISE_One_NoSignal"
        return "IMPRECISE"
    elif filter_val == "PRECISE":
        merged_hap = info_dict.get("MERGED_HAP")
        if merged_hap == "PRECISE_CONF":
            return "PRECISE_CONF"
        elif merged_hap == "PRECISE_Mid_CONF":
            return "PRECISE_Mid_CONF"
        elif merged_hap == "PRECISE_Low_CONF":
            return "PRECISE_Low_CONF"
        return "PRECISE"
    return "UNKNOWN"

def parse_vcf_step7(filepath, caller_name):
    """
    Parses a VCF file and extracts relevant information. (Step 7)
    """
    variants_data = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    print(f"(Step 7) Skipping malformed line: {line.strip()}", file=sys.stderr)
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                filter_val = fields[6]
                info_str = fields[7]
                info_dict = parse_info_field_step7(info_str)

                original_locus = info_dict.get("ORIGINAL_LOCUS")
                if not original_locus:
                    continue

                category = determine_category_step7(filter_val, info_dict)

                variants_data.append({
                    "ORIGINAL_LOCUS": original_locus,
                    f"{caller_name}_category": category,
                    f"{caller_name}_CHROM": chrom,
                    f"{caller_name}_POS": pos,
                    f"{caller_name}_FILTER": filter_val,
                })
    except FileNotFoundError:
        print(f"Error (Step 7): Input VCF not found: {filepath}", file=sys.stderr)
        return pd.DataFrame()
    except Exception as e:
        print(f"Error (Step 7) parsing VCF {filepath}: {e}", file=sys.stderr)
        return pd.DataFrame()

    return pd.DataFrame(variants_data)

def resolve_duplicates_step7(variants_df, caller_name):
    """
    Resolves multiple calls for the same ORIGINAL_LOCUS. (Step 7)
    """
    if variants_df.empty:
        return pd.DataFrame(columns=['ORIGINAL_LOCUS', f'{caller_name}_category', f'{caller_name}_CHROM', f'{caller_name}_POS', f'{caller_name}_FILTER'])

    priority_order = {
        "PRECISE_CONF": 1,
        "PRECISE_Mid_CONF": 2,
        "PRECISE_Low_CONF": 3,
        "PRECISE": 4,
        "IMPRECISE_Two_Hap": 5,
        "IMPRECISE_One_NoSignal": 6,
        "IMPRECISE": 7,
        "NOSIGNAL": 8,
        "UNKNOWN": 9
    }

    variants_df['priority'] = variants_df[f'{caller_name}_category'].map(priority_order).fillna(99)
    resolved_df = variants_df.sort_values(by=['ORIGINAL_LOCUS', 'priority']).drop_duplicates(subset=['ORIGINAL_LOCUS'], keep='first')
    resolved_df = resolved_df.drop(columns=['priority'])
    return resolved_df

def compare_vcfs_advanced_step7(vcf_file1_path, vcf_file2_path, caller1_name, caller2_name, output_log_file, output_summary_file):
    """
    Core logic for Step 7.
    """
    print(f"(Step 7) Processing {caller1_name} VCF: {vcf_file1_path}")
    df1_raw = parse_vcf_step7(vcf_file1_path, caller1_name)
    print(f"(Step 7) Found {len(df1_raw)} raw variants in {vcf_file1_path}")

    print(f"(Step 7) Processing {caller2_name} VCF: {vcf_file2_path}")
    df2_raw = parse_vcf_step7(vcf_file2_path, caller2_name)
    print(f"(Step 7) Found {len(df2_raw)} raw variants in {vcf_file2_path}")

    df1_resolved = resolve_duplicates_step7(df1_raw, caller1_name)
    print(f"(Step 7) Resolved {len(df1_resolved)} unique loci for {caller1_name}")
    df2_resolved = resolve_duplicates_step7(df2_raw, caller2_name)
    print(f"(Step 7) Resolved {len(df2_resolved)} unique loci for {caller2_name}")

    merged_df = pd.merge(df1_resolved, df2_resolved, on="ORIGINAL_LOCUS", how="outer", suffixes=(f"_{caller1_name}", f"_{caller2_name}"))

    log_data = []
    comparison_data_for_table = []

    defined_categories = [
        "IMPRECISE", "IMPRECISE_One_NoSignal", "IMPRECISE_Two_Hap",
        "PRECISE_CONF", "PRECISE_Mid_CONF", "PRECISE_Low_CONF",
        "NOSIGNAL"
    ]
    
    row_categories = sorted(list(set(defined_categories)))
    col_categories = sorted(list(set(defined_categories))) + ["NOT_FOUND"]


    for _, row in merged_df.iterrows():
        locus = row['ORIGINAL_LOCUS']
        cat1 = row.get(f'{caller1_name}_category', 'NOT_FOUND')
        cat2 = row.get(f'{caller2_name}_category', 'NOT_FOUND')

        details1 = f"CHROM:{row.get(f'{caller1_name}_CHROM', 'N/A')},POS:{row.get(f'{caller1_name}_POS', 'N/A')},FILTER:{row.get(f'{caller1_name}_FILTER', 'N/A')}"
        details2 = f"CHROM:{row.get(f'{caller2_name}_CHROM', 'N/A')},POS:{row.get(f'{caller2_name}_POS', 'N/A')},FILTER:{row.get(f'{caller2_name}_FILTER', 'N/A')}"
        log_data.append(f"{locus}\t{cat1}\t{cat2}\t{details1}\t{details2}\n")

        final_cat1_for_table = cat1
        final_cat2_for_table = cat2

        if cat1 == "NOT_FOUND": 
            continue
        
        if cat1 == "NOSIGNAL":
            if cat2 == "NOSIGNAL":
                final_cat2_for_table = "NOSIGNAL"
            elif cat2 == "NOT_FOUND":
                final_cat2_for_table = "NOT_FOUND"

        if final_cat1_for_table != "NOT_FOUND":
            comparison_data_for_table.append({
                f"{caller1_name}_Category": final_cat1_for_table,
                f"{caller2_name}_Category": final_cat2_for_table
            })

    try:
        with open(output_log_file, 'w') as log_f:
            log_f.write(f"ORIGINAL_LOCUS\t{caller1_name}_Category\t{caller2_name}_Category\t{caller1_name}_Details\t{caller2_name}_Details\n")
            for entry in log_data:
                log_f.write(entry)
    except IOError as e:
        print(f"Error (Step 7) writing log file {output_log_file}: {e}", file=sys.stderr)


    if not comparison_data_for_table:
        print("(Step 7) No data to create a summary table.")
        with open(output_summary_file, 'w') as f_summary:
            f_summary.write("No comparable data found based on the specified criteria.\n")
        return

    summary_df = pd.DataFrame(comparison_data_for_table)

    summary_df = summary_df[summary_df[f'{caller1_name}_Category'].isin(row_categories)]
    summary_df = summary_df[summary_df[f'{caller2_name}_Category'].isin(col_categories)]
    
    if summary_df.empty:
        print("(Step 7) No data to create a summary table after filtering for defined categories.")
        with open(output_summary_file, 'w') as f_summary:
            f_summary.write("No comparable data found after filtering for defined categories.\n")
        return

    summary_table = pd.crosstab(
        pd.Categorical(summary_df[f'{caller1_name}_Category'], categories=row_categories, ordered=True),
        pd.Categorical(summary_df[f'{caller2_name}_Category'], categories=col_categories, ordered=True),
        dropna=False
    )

    summary_table.index.name = f"{caller1_name}_Category"
    summary_table.columns.name = f"{caller2_name}_Category"

    print("\n(Step 7) Comparison Summary Table:")
    print(summary_table)

    try:
        summary_table.to_csv(output_summary_file, sep='\t')
        print(f"\n(Step 7) Summary table saved to: {output_summary_file}")
    except IOError as e:
        print(f"Error (Step 7) writing summary file {output_summary_file}: {e}", file=sys.stderr)


def run_step7(vcf_file1, vcf_file2, name1, name2, output_log, output_summary):
    """
    Runner function for Step 7.
    """
    try:
        compare_vcfs_advanced_step7(vcf_file1, vcf_file2, name1, name2, output_log, output_summary)
        print(f"(Step 7) Comparison complete.")
        print(f"(Step 7) Detailed log saved to: {output_log}")
    except Exception as e:
        print(f"Error in Step 7 (table_matrix) for {vcf_file1} / {vcf_file2}: {e}", file=sys.stderr)
        sys.exit(1)


# ############################################################################
# --- MAIN PIPELINE ORCHESTRATOR ---
# ############################################################################

def main_pipeline():
    """
    Main orchestrator function for the integrated pipeline.
    Parses all arguments and calls steps in sequence.
    """
    parser = argparse.ArgumentParser(
        description="Integrated pipeline for inversion breakpoint refinement from Strand-seq BAMs.",
        formatter_class=argparse.RawTextHelpFormatter
    )
    
    # --- Step 1 Inputs ---
    g1 = parser.add_argument_group("Step 1: BAM to VCF Inputs")
    g1.add_argument("--bed", required=True, help="Input BED file with inversion locations (for Step 1 & 3).")
    g1.add_argument("--p_bam_mm2", required=True, help="Paternal BAM file (Aligner 1: minimap2).")
    g1.add_argument("--m_bam_mm2", required=True, help="Maternal BAM file (Aligner 1: minimap2).")
    g1.add_argument("--p_bam_vac", required=True, help="Paternal BAM file (Aligner 2: VACmap).")
    g1.add_argument("--m_bam_vac", required=True, help="Maternal BAM file (Aligner 2: VACmap).")
    g1.add_argument("-t", "--threads", type=int, default=8, help="Number of CPU threads (for Step 1).")
    g1.add_argument("--debug", action="store_true", help="Enable debug mode (for Step 1).")

    # --- Step 4 Inputs ---
    g4 = parser.add_argument_group("Step 4: Annotation Inputs")
    g4.add_argument("--bed_5kb", required=True, help="Path to repeats_ge5kb.bed (for Step 4).")
    g4.add_argument("--bed_10kb", required=True, help="Path to repeats_ge10kb.bed (for Step 4).")
    g4.add_argument("--bed_centromeres", required=True, help="Path to centromeres_hg38.bed (for Step 4).")

    # --- Output Naming ---
    g_out = parser.add_argument_group("Output Configuration")
    g_out.add_argument("-o", "--output_prefix", default="pipeline_run", 
                         help="Prefix for all intermediate and final files (e.g., /path/to/my_run).")

    args = parser.parse_args()
    
    # Define file paths based on prefix
    # Ensure the directory for the prefix exists
    work_dir = os.path.abspath(os.path.dirname(args.output_prefix))
    prefix = os.path.basename(args.output_prefix)
    
    if not os.path.exists(work_dir):
        try:
            os.makedirs(work_dir, exist_ok=True)
            print(f"Created output directory: {work_dir}")
        except OSError as e:
            print(f"FATAL: Could not create output directory {work_dir}: {e}", file=sys.stderr)
            sys.exit(1)
    
    # Helper to create paths
    def p(step, name, ext):
        return os.path.join(work_dir, f"{prefix}.step{step}.{name}.{ext}")

    # --- File Path Definitions ---
    s1_log_mm2 = p(1, "mm2_bam2vcf", "log")
    s1_out_mm2 = p(1, "mm2", "vcf")
    s1_log_vac = p(1, "vacmap_bam2vcf", "log")
    s1_out_vac = p(1, "vacmap", "vcf")
    
    s2_out_mm2 = p(2, "mm2.merged", "vcf")
    s2_out_vac = p(2, "vacmap.merged", "vcf")
    
    s3_log     = p(3, "rm_redundant", "log")
    s3_out_mm2 = p(3, "mm2.filtered", "vcf")
    s3_out_vac = p(3, "vacmap.filtered", "vcf")
    
    s4_out_mm2 = p(4, "mm2.annotated", "vcf")
    s4_out_vac = p(4, "vacmap.annotated", "vcf")
    
    s5_out_final = p(5, "final_merged", "vcf")
    # s5 also creates .dropped.vcf and .log based on s5_out_final
    
    s6_out_stats = p(6, "final_stats", "txt")
    
    s7_out_log     = p(7, "comparison", "log.tsv")
    s7_out_summary = p(7, "comparison_summary", "tsv")

    print(f"--- Pipeline Starting ---")
    print(f"All outputs will be prefixed with: {args.output_prefix}")

    # --- RUN STEP 1 (mm2) ---
    print(f"\n[Step 1/7] Running bam2vcf for minimap2...")
    step1_args_mm2 = {
        "bed": args.bed,
        "paternal_bam": args.p_bam_mm2,
        "maternal_bam": args.m_bam_mm2,
        "output_vcf": s1_out_mm2,
        "debug": args.debug,
        "threads": args.threads,
        "log_file": s1_log_mm2
    }
    run_step1(step1_args_mm2)
    print(f"Finished Step 1 (minimap2). Output: {s1_out_mm2}")

    # --- RUN STEP 1 (vac) ---
    print(f"\n[Step 1/7] Running bam2vcf for VACmap...")
    step1_args_vac = {
        "bed": args.bed,
        "paternal_bam": args.p_bam_vac,
        "maternal_bam": args.m_bam_vac,
        "output_vcf": s1_out_vac,
        "debug": args.debug,
        "threads": args.threads,
        "log_file": s1_log_vac
    }
    run_step1(step1_args_vac)
    print(f"Finished Step 1 (VACmap). Output: {s1_out_vac}")

    # --- RUN STEP 2 (mm2) ---
    print(f"\n[Step 2/7] Merging haplotypes for minimap2...")
    run_step2(s1_out_mm2, s2_out_mm2)
    print(f"Finished Step 2 (minimap2). Output: {s2_out_mm2}")

    # --- RUN STEP 2 (vac) ---
    print(f"\n[Step 2/7] Merging haplotypes for VACmap...")
    run_step2(s1_out_vac, s2_out_vac)
    print(f"Finished Step 2 (VACmap). Output: {s2_out_vac}")

    # --- RUN STEP 3 ---
    print(f"\n[Step 3/7] Removing redundant inversions...")
    run_step3(args.bed, s2_out_mm2, s2_out_vac, s3_log, s3_out_mm2, s3_out_vac)
    print(f"Finished Step 3. Outputs: {s3_out_mm2}, {s3_out_vac}")

    # --- RUN STEP 4 (mm2) ---
    print(f"\n[Step 4/7] Annotating minimap2 VCF...")
    run_step4(s3_out_mm2, s4_out_mm2, args.bed_5kb, args.bed_10kb, args.bed_centromeres)
    print(f"Finished Step 4 (minimap2). Output: {s4_out_mm2}")

    # --- RUN STEP 4 (vac) ---
    print(f"\n[Step 4/7] Annotating VACmap VCF...")
    run_step4(s3_out_vac, s4_out_vac, args.bed_5kb, args.bed_10kb, args.bed_centromeres)
    print(f"Finished Step 4 (VACmap). Output: {s4_out_vac}")

    # --- RUN STEP 5 ---
    print(f"\n[Step 5/7] Merging aligner VCFs...")
    run_step5(s4_out_mm2, s4_out_vac, s5_out_final)
    print(f"Finished Step 5. Output: {s5_out_final}")

    # --- RUN STEP 6 ---
    print(f"\n[Step 6/7] Generating stats for final VCF...")
    run_step6(s5_out_final, s6_out_stats)
    print(f"Finished Step 6. Output: {s6_out_stats}")

    # --- RUN STEP 7 ---
    print(f"\n[Step 7/7] Comparing annotated aligner VCFs...")
    run_step7(s4_out_mm2, s4_out_vac, "minimap2", "VACmap", s7_out_log, s7_out_summary)
    print(f"Finished Step 7. Output: {s7_out_summary}")

    print(f"\n--- Pipeline Finished Successfully ---")
    print(f"Final merged VCF: {s5_out_final}")
    print(f"Final stats: {s6_out_stats}")
    print(f"Aligner comparison: {s7_out_summary}")


if __name__ == "__main__":
    # For Step 1's multiprocessing
    multiprocessing.freeze_support()
    main_pipeline()
