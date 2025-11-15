#!/usr/bin/env bash
set -euo pipefail

#############################
# CONFIGURE THESE TWO PATHS #
#############################

# Folder that contains the *calling results*
CALL_BASE="/users/u250191/ryan_scratch_ln/benchmark_inv/round2_analysis_251111/different_coverage/caller"

# Folder that contains the *destination* structure
DEST_BASE="/users/u250191/ryan_scratch_ln/benchmark_inv/round2_analysis_251111/script_for_github/data_zenodo/inv_call_5samples"

#############################
# CONSTANTS / ARRAYS        #
#############################

coverages=(5x 10x 20x 30x)
samples=(hg002 hg00733 hg02818 hg03486 na19240)

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

    echo "=== Updating Manta for ${cov}/${sample} ==="

    ########################
    # Manta Only           #
    ########################

    # Source changed from inv_only.vcf to inv_with_infor.vcf
    manta_src="${src_root}/manta/inv_with_infor.vcf"
    
    # Destination remains the same standard name
    manta_dest="${dest_root}/manta_inv.vcf"
    
    safe_cp "$manta_src" "$manta_dest"

  done   # sample
done     # cov
