#!/usr/bin/env bash
set -euo pipefail
root="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$root/Data/raw/1000g"
echo "Download the 1000 Genomes phase-3 EUR PLINK reference panel to Data/raw/1000g."
