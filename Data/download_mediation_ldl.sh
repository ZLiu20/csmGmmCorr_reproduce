#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$root/Data/raw/mediation_ldl"
cat <<'EOF'
Place QTD000266.all.tsv.gz, gencode.v26.GRCh38.genes.gtf,
GTEx_v8_lookup_table.txt.gz, mashr_liver.db and mashr_liver.txt.gz in
Data/raw/mediation_ldl. These public resources remain subject to their source
repositories' terms and are therefore not silently mirrored here.
EOF
