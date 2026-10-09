#!/bin/bash
# Authors' step, done once: store the values of the paper in reference/ (used by tests/verify.R).
# Run from the repository root of the working copy that holds the results of the paper:
#   bash tests/snapshot_reference.sh            make all figures and tables, then copy them
#   bash tests/snapshot_reference.sh --link     first link the results of the original HPC runs
#                                               into the folder layout (author's folders)
# Writes reference/ and reference.tar.gz. Files are copied, never linked.
set -u
cd "$(dirname "$0")/.."
module load R/4.5.1 2>/dev/null; mkdir -p logs
if [ "${1:-}" = "--link" ]; then
  c=$(git log --diff-filter=D --format=%H -1 -- tests/link_existing_results.sh)
  git show "$c^:tests/link_existing_results.sh" > /tmp/link_existing_results_$$.sh
  bash /tmp/link_existing_results_$$.sh; rm -f /tmp/link_existing_results_$$.sh
fi
[ -s logs/versions.txt ] || Rscript setup.R > logs/setup.log 2>&1
VALIDATION=0 bash make_all.sh
rm -rf reference; mkdir -p reference
copy() { for f in "$@"; do [ -f "$f" ] || continue; mkdir -p "reference/$(dirname "$f")"; cp -L "$f" "reference/$f"; done; }
shopt -s nullglob
copy Fig*/*.tex Tab*/*.tex SuppTab*/*.tex Fig*/*.csv Tab*/*.csv SuppTab*/*.csv Fig1/Fig1_numbers.txt Fig4/Fig4_numbers.txt
copy Fig2/output/*.csv Fig3/output/*.csv Fig4/output/*.csv SuppTab4/output/*.csv SuppTab4/output/*.tsv
copy Fig4/output/*.tsv Fig4/output/mediation_ldl/compare/*.tsv
# the 13 application fits and the 23 mediation fits of RD_JOBS and MEDLDL_JOBS only
for d in $(Rscript -e 'here::i_am("setup.R"); source(here::here("SupportingCode","load_all.R")); cat(sub(paste0(here::here(), "/"), "", c(rd_fit_dirs(), medldl_fit_dirs()), fixed = TRUE))' 2>/dev/null | tail -1); do
  copy $d/summary.tsv $d/sigma.tsv $d/all_t.tsv
done
cp logs/versions.txt reference/versions.txt
cat > reference/README.md <<EOF
Values of the paper: the figures and tables made on $(date +%F) by make_all.sh from the results of
the original runs (tests/snapshot_reference.sh). versions.txt lists the R and package versions
installed when the reference was made. tests/verify.R compares a new run with these files.
EOF
find reference -type f | sort > reference/files.txt
echo "$(wc -l < reference/files.txt) files in reference/"
grep -E "FAILED" logs/make/status.txt
Rscript tests/verify.R > /dev/null; head -12 logs/verify_report.txt     # self-check: all IDENTICAL
tar czf reference.tar.gz reference && ls -lh reference.tar.gz
