#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$root/Data/raw/glgc/pleiotropy/GLGC_trans" "$root/Data/raw/glgc/replication/GLGC_LDL_sex_ancestry" "$root/Data/raw/ukb"
cat <<'EOF'
Download the GLGC 2021 trans-ancestry and ancestry-by-sex files named in
Data/glgc_pleiotropy_manifest.tsv and Data/glgc_replication_manifest.tsv.
UK Biobank Neale Lab files for traits 21001, 4080, 30870, 30760 and 30750
belong in Data/raw/ukb. Source URLs and access conditions may change, so this
script creates the required layout without bypassing provider terms.
EOF
