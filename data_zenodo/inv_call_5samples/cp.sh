#!/usr/bin/env bash
set -euo pipefail

#############################
# CONFIGURE THESE TWO PATHS #
#############################

# Folder that contains the *calling results*:
#   5x/hg002/cuteSV2/, 5x/hg002/gridss2/, 5x/hg002/manta/, ...
CALL_BASE="/users/u250191/ryan_scratch_ln/benchmark_inv/round2_analysis_251111/different_coverage/caller"

# Folder that contains the *destination* structure:
#   5x/hg002/, 10x/hg00733/, etc. (where you want to copy/rename VCFs to)
DEST_BASE="/users/u250191/ryan_scratch_ln/benchmark_inv/round2_analysis_251111/script_for_github/data_zenodo/inv_call_5samples"

#############################
# CONSTANTS / ARRAYS       #
#############################

coverages=(5x 10x 20x 30x)
samples=(hg002 hg00733 hg02818 hg03486 na19240)

# Long-read callers that have per-combination subfolders
lr_callers=(cuteSV2 severus sniffles2 svim svision-pro)
combs=(ont_minimap2 ont_vacmap pb_minimap2 pb_vacmap)

#############################
# HELPER: safe_cp           #
#############################
safe_cp() {
    local src="$1"
    local dest="$2"

    if [[ -f "$src" ]]; then
        echo "Copying: $src  ->  $dest"
        cp "$src" "$dest"
    else
        echo "WARNING: missing file: $src" >&2
    fi
}

#############################
# MAIN LOOP                 #
#############################

for cov in "${coverages[@]}"; do
  for sample in "${samples[@]}"; do
    src_root="${CALL_BASE}/${cov}/${sample}"
    dest_root="${DEST_BASE}/${cov}/${sample}"

    # Skip if caller root doesn't exist
    if [[ ! -d "$src_root" ]]; then
        echo "WARNING: missing source dir: $src_root" >&2
        continue
    fi

    # Make sure destination exists
    mkdir -p "$dest_root"

    echo "=== Processing ${cov}/${sample} ==="

    ########################
    # Short-read callers   #
    ########################

    # gridss2: simple_inv.vcf -> gridss2_inv.vcf
    gridss_src="${src_root}/gridss2/simple_inv.vcf"
    gridss_dest="${dest_root}/gridss2_inv.vcf"
    safe_cp "$gridss_src" "$gridss_dest"

    # manta: inv_only.vcf -> manta_inv.vcf
    manta_src="${src_root}/manta/inv_only.vcf"
    manta_dest="${dest_root}/manta_inv.vcf"
    safe_cp "$manta_src" "$manta_dest"

    ########################
    # Long-read callers    #
    ########################

    for comb in "${combs[@]}"; do
      comb_dir_cutesv2="${src_root}/cuteSV2/${comb}"
      comb_dir_severus="${src_root}/severus/${comb}"
      comb_dir_sniffles="${src_root}/sniffles2/${comb}"
      comb_dir_svim="${src_root}/svim/${comb}"
      comb_dir_svision="${src_root}/svision-pro/${comb}"

      # ---- cuteSV2 ----
      # Pattern:
      #   - hg002 can have: inv_cutesv2_989518.vcf
      #   - others:         inv_cutesv2_ont_minimap2.vcf, etc.
      # We just grab the first inv_cutesv2*.vcf we see in the comb folder.
      if [[ -d "$comb_dir_cutesv2" ]]; then
          cutesv_src=$(ls "${comb_dir_cutesv2}"/inv_cutesv2*.vcf 2>/dev/null | head -n1 || true)
          if [[ -n "${cutesv_src:-}" ]]; then
              cutesv_dest="${dest_root}/inv_cutesv2_${comb}.vcf"
              safe_cp "$cutesv_src" "$cutesv_dest"
          else
              echo "WARNING: no inv_cutesv2*.vcf in ${comb_dir_cutesv2}" >&2
          fi
      fi

      # ---- severus ----
      # Pattern: add_end_inv_severus.vcf
      if [[ -d "$comb_dir_severus" ]]; then
          severus_src="${comb_dir_severus}/add_end_inv_severus.vcf"
          severus_dest="${dest_root}/inv_severus_${comb}.vcf"
          safe_cp "$severus_src" "$severus_dest"
      fi

      # ---- sniffles2 ----
      # Good file:    inv_sniffle2_989527.vcf
      # Bad/ignore:   inv_sniffle2_989886_woheader.vcf
      if [[ -d "$comb_dir_sniffles" ]]; then
          sniff_src=$(ls "${comb_dir_sniffles}"/inv_sniffle2*.vcf 2>/dev/null | grep -v "woheader" | head -n1 || true)
          if [[ -n "${sniff_src:-}" ]]; then
              sniff_dest="${dest_root}/inv_sniffles2_${comb}.vcf"
              safe_cp "$sniff_src" "$sniff_dest"
          else
              echo "WARNING: no usable inv_sniffle2*.vcf (without 'woheader') in ${comb_dir_sniffles}" >&2
          fi
      fi

      # ---- svim ----
      # Pattern: filtered_qual_10_inv.vcf
      if [[ -d "$comb_dir_svim" ]]; then
          svim_src="${comb_dir_svim}/filtered_qual_10_inv.vcf"
          svim_dest="${dest_root}/inv_svim_${comb}.vcf"
          safe_cp "$svim_src" "$svim_dest"
      fi

      # ---- svision-pro ----
      # Pattern: inv_only_add_simplesv.sorted.vcf
      # Naming:  inv_svisionpro_<comb>.vcf
      if [[ -d "$comb_dir_svision" ]]; then
          svision_src="${comb_dir_svision}/inv_only_add_simplesv.sorted.vcf"
          svision_dest="${dest_root}/inv_svisionpro_${comb}.vcf"
          safe_cp "$svision_src" "$svision_dest"
      fi

    done # comb
  done   # sample
done     # cov

