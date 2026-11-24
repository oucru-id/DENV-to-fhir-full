#!/usr/bin/env bash
set -euo pipefail

# Pre-fetch all 4 Dengue Nextclade datasets for offline/reproducible DENVtoFHIR runs.
# Usage: tools/nextclade/fetch_datasets.sh [output_dir] [optional_tag]

OUT_DIR="${1:-db/nextclade}"
TAG="${2:-}"
TAG_ARG=()
if [ -n "$TAG" ]; then
    TAG_ARG=(--tag "$TAG")
fi

mkdir -p "$OUT_DIR"
for sero in denv1 denv2 denv3 denv4; do
    echo "Fetching community/v-gen-lab/dengue/${sero} -> ${OUT_DIR}/${sero}"
    nextclade dataset get \
        --name "community/v-gen-lab/dengue/${sero}" \
        "${TAG_ARG[@]}" \
        --output-dir "${OUT_DIR}/${sero}"
done

echo "Done. Pass --nextclade_dir ${OUT_DIR} to DENVtoFHIR."
