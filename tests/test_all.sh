#!/bin/bash
# Test of the repository on the cluster (authors' check before release). Run from the
# repository root:
#   bash tests/test_all.sh local     syntax check, setup, link the original results, and run every
#                                    plot_* / make_* script (login node; about 30-60 minutes)
#   bash tests/test_all.sh submit    LSF smoke jobs: task 1 of every simulation array with
#                                    SIM_N = 2000 and one replicate; real-data fit on chromosome 22;
#                                    the weight-floor study (full run); R port of the GLGC preparation
#   bash tests/test_all.sh report    collect the logs, PDFs and tables into tests/report.tar.gz
set -u
mode=${1:-local}
module load R/4.5.1 2>/dev/null
case $mode in
local)
  Rscript tests/check_syntax.R || exit 1
  Rscript setup.R ${CSMGMMCORR_SRC:-} 2>&1 | tee tests/setup.log | tail -15
  bash tests/link_existing_results.sh | tee tests/link.log
  mkdir -p tests/logs
  for s in Fig1/plot_fig1.R Fig2/plot_fig2.R SuppFig1/plot_suppfig1_12.R Tab1/make_tab1.R Fig3/plot_fig3.R Tab2/make_tab2.R \
           SuppTab2/make_supptab2.R SuppTab3/make_supptab3.R SuppTab5/make_supptab5.R SuppTab6/make_supptab6_8.R \
           SuppTab9/make_supptab9.R SuppTab10/make_supptab10_11.R SuppTab12/make_supptab12_15.R SuppTab16/make_supptab16_18.R \
           Fig4/make_tab3.R Fig4/plot_fig4.R Fig4/plot_fig5.R Fig4/holdout_validation.R SuppTab19/make_supptab19_23.R; do
    log=tests/logs/$(basename $s .R).log
    if Rscript $s > $log 2>&1; then echo "OK      $s"; else echo "FAILED  $s   (see $log)"; tail -5 $log; fi
  done
  Rscript -e 'here::i_am("setup.R"); source(here::here("SupportingCode","load_all.R")); print(sim_progress(), row.names = FALSE)' > tests/logs/progress.log 2>&1
  ;;
submit)
  mkdir -p tests/errDir
  n=$(Rscript -e 'here::i_am("setup.R"); source(here::here("SupportingCode","load_all.R")); cat(length(sim_arrays()))' 2>/dev/null | tail -1)
  bsub -J "smoke[1-$n]" -q medium -n 4 -W 4:00 -M 8 -R "span[hosts=1] rusage[mem=8]" -o "tests/errDir/smoke_%I.out" -e "tests/errDir/smoke_%I.err" <<JOB
module load R/4.5.1
export SIM_N=2000 MAX_ITER=200 SIM_MAX_REPS=1 OMP_NUM_THREADS=1
Rscript -e 'here::i_am("setup.R"); source(here::here("SupportingCode","load_all.R")); a <- names(sim_arrays())[as.integer(Sys.getenv("LSB_JOBINDEX"))]; message("array ", a); run_sim_task(a, 1L)'
JOB
  bsub -J "smoke_rd" -q medium -n 8 -W 4:00 -M 32 -R "span[hosts=1] rusage[mem=32]" -o tests/errDir/smoke_rd.out -e tests/errDir/smoke_rd.err <<JOB
module load R/4.5.1
export SMOKE=1
Rscript Fig4/realdata_fit.R 1
JOB
  bsub < SuppTab4/run_SuppTab4.lsf
  RAW=${RAW_GLGC:-$(dirname $(dirname $(find /rsrch8/home/biostatistics/zliu20 -maxdepth 8 -name "LDL_INV_EUR_HRC_1KGP3_others_MALE.meta.singlevar.results.gz" 2>/dev/null | head -1)))}
  echo "raw GLGC folder: $RAW"
  bsub -J "glgc_prep" -q medium -n 1 -W 8:00 -M 64 -R "span[hosts=1] rusage[mem=64]" -o tests/errDir/glgc_prep.out -e tests/errDir/glgc_prep.err <<JOB
module load R/4.5.1
Rscript tests/compare_glgc_prep.R "$RAW" "$(readlink -f Data/prepared/glgc_replication)"
JOB
  ;;
report)
  Rscript SuppTab4/make_supptab4.R > tests/logs/make_supptab4.log 2>&1
  for f in tests/errDir/smoke_*.err; do grep -l -i "error" $f; done > tests/logs/smoke_errors.txt 2>/dev/null
  ls Fig2/output_smoke/*/* Fig3/output_smoke/*/* Tab2/output_smoke/*/* SuppTab*/output_smoke/*/* 2>/dev/null | wc -l > tests/logs/smoke_count.txt
  tar czf tests/report.tar.gz tests/*.log tests/logs tests/errDir/smoke_rd.* tests/errDir/glgc_prep.* \
      Fig1/*.pdf Fig1/*.txt Fig2/*.pdf Fig3/*.pdf Fig4/*.pdf Fig4/*.tex Fig4/*.txt Fig4/output/holdout_validation.csv SuppFig1/*.pdf \
      Tab1/*.tex Tab2/*.tex SuppTab*/*.tex SuppTab*/*.csv SuppTab4/output 2>/dev/null
  ls -lh tests/report.tar.gz
  ;;
esac
