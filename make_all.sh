#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
module load R/4.5.1 2>/dev/null || true
mkdir -p logs/make
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
rm -f */output/all_results_*.csv.gz

scripts=(
  Fig1/plot_fig1.R Fig2/plot_fig2.R SuppFig1/plot_suppfig1_12.R
  Tab1/make_tab1.R Fig3/plot_fig3.R Tab2/make_tab2.R
  SuppTab2/make_supptab2.R SuppTab3/make_supptab3.R
  SuppTab4/make_supptab4.R SuppTab5/make_supptab5.R
  SuppTab6/make_supptab6_8.R SuppTab9/make_supptab9.R
  SuppTab10/make_supptab10_11.R SuppTab12/make_supptab12_15.R
  SuppTab16/make_supptab16_18.R Fig4/make_tab3.R Fig4/plot_fig4.R
  Fig4/plot_fig5.R Fig4/holdout_validation.R SuppTab19/make_supptab19_23.R
  Fig4/summarize_mediation_ldl.R Fig4/diagnose_mediation_ldl.R
  SuppTab24/make_supptab24_26.R
)
if [[ "${VALIDATION:-0}" == 1 ]]; then
  scripts+=(Validation/validate_csmgmm.R Validation/validate2_csmgmm.R Validation/diag_incon.R)
fi

: > logs/make/status.txt
failed=0
for script in "${scripts[@]}"; do
  log="logs/make/$(basename "${script%.R}").log"
  started=$(date +%s)
  if Rscript "$script" >"$log" 2>&1; then state=OK; else state=FAILED; failed=$((failed+1)); fi
  printf '%-36s %-7s %6ss %s\n' "$script" "$state" "$(( $(date +%s) - started ))" "$log" | tee -a logs/make/status.txt
  [[ "$state" == FAILED ]] && tail -n 5 "$log"
done
printf 'make_all: %d script(s) failed\n' "$failed" | tee -a logs/make/status.txt
exit 0
