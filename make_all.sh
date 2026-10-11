#!/bin/bash
# Make every figure and table from the finished computations (run from the repository root).
# Called by the final job of run_all.sh; can also be run by hand:   bash make_all.sh
# One log per script in logs/make/; a script that fails does not stop the others.
# VALIDATION=0 skips the package checks of Validation/ (they install and run Sun et al.'s csmGmm).
set -u
cd "$(dirname "$0")"
module load R/4.5.1 2>/dev/null
mkdir -p logs/make; : > logs/make/status.txt
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
# combined simulation tables are cached by read_sim(); rebuild them from the replicate files
rm -f */output/all_results_*.csv.gz
SCRIPTS="Fig1/plot_fig1.R
Fig2/plot_fig2.R
SuppFig1/plot_suppfig1_12.R
Tab1/make_tab1.R
Fig3/plot_fig3.R
Tab2/make_tab2.R
SuppTab2/make_supptab2.R
SuppTab3/make_supptab3.R
SuppTab4/make_supptab4.R
SuppTab5/make_supptab5.R
SuppTab6/make_supptab6_8.R
SuppTab9/make_supptab9.R
SuppTab10/make_supptab10_11.R
SuppTab12/make_supptab12_15.R
SuppTab16/make_supptab16_18.R
Fig4/make_tab3.R
Fig4/plot_fig4.R
Fig4/plot_fig5.R
Fig4/holdout_validation.R
SuppTab19/make_supptab19_23.R
Fig4/summarize_mediation_ldl.R
Fig4/diagnose_mediation_ldl.R
SuppTab24/make_supptab24_26.R"
[ "${VALIDATION:-1}" = "1" ] && SCRIPTS="$SCRIPTS
Validation/validate_csmgmm.R
Validation/validate2_csmgmm.R
Validation/diag_incon.R"
nf=0
for s in $SCRIPTS; do
  log=logs/make/$(basename "$s" .R).log
  t0=$(date +%s)
  if Rscript "$s" > "$log" 2>&1; then st=OK; else st=FAILED; nf=$((nf + 1)); fi
  printf "%-7s %-36s %5ds   %s\n" "$st" "$s" $(( $(date +%s) - t0 )) "$log" | tee -a logs/make/status.txt
  [ $st = FAILED ] && tail -n 5 "$log" | sed 's/^/        /'
done
echo "$(date '+%F %T') make_all: $nf script(s) failed" | tee -a logs/make/status.txt
exit 0
