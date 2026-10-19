#!/bin/bash
# Reproduce every figure and table of the paper and the Supplementary Material in one run.
# Run from the repository root of a fresh clone, on a login node with internet access:
#
#   nohup bash run_all.sh > logs/run_all.log 2>&1 &
#
# Steps (each is skipped if already done, so the command can simply be repeated to resume):
#   1. setup.R: packages, compiled kernels, logs/versions.txt
#   2. submit every simulation array (Fig2, Fig3, Tab2, SuppTab3, SuppTab4, SuppTab9, SuppTab12);
#      they need no data and start at once
#   3. download the summary statistics (GLGC, UK Biobank; about 25 GB) while the simulations run
#   4. submit the data preparation and the application fits with LSF dependencies:
#        prepare[1-3]  -> realdata[1-13]
#        medscc[1-11], medscc_cap[1-4]   (lung SCC mediation; only if the non-public data have
#                                         been prepared, see Data/README.md)
#   5. submit one final job that waits for all of the above and runs make_all.sh
# Dependencies use ended() rather than done(): a failed step does not leave the later jobs
# pending for ever; they fail too, and logs/make/status.txt lists the scripts that failed. Steps whose
# results are complete (all replicates of an array, DONE files of the fits) are not submitted
# again, so after fixing a failure the same command reruns only what is missing.
# The figures and tables are written to the folders Fig*/, Tab*/, SuppFig*/ and SuppTab*/;
# logs/make/status.txt lists every script with OK or FAILED.
#
# Other modes:
#   bash run_all.sh status      progress of the simulation arrays and state of the LSF jobs
#   bash run_all.sh final       submit only the final job (e.g. after resubmitting failed tasks)
#   bash run_all.sh make        run make_all.sh here (interactive)
#   DRYRUN=1 bash run_all.sh    show what would be submitted (no jobs, no downloads)
# Environment: QUEUE (default medium) for the jobs submitted here; SKIP_DOWNLOAD=1 if Data/raw is
# already complete (for example a symbolic link to an earlier download); CSMGMMCORR_SRC=<path to
# csmGmmCorr_1.3.0.tar.gz or its source folder> to install the package from a local copy.
set -u
cd "$(dirname "$0")"
mkdir -p logs
QUEUE=${QUEUE:-medium}
if [ "${DRYRUN:-0}" = "1" ]; then bsub() { cat > /dev/null; echo "[dry run] bsub $*"; }; fi
module load R/4.5.1 2>/dev/null
mode=${1:-all}
stamp() { echo "[$(date '+%F %T')] $*"; }

jobname() { sed -n 's/^#BSUB -J \([A-Za-z0-9_]*\).*/\1/p' "$1" | head -1; }
errdir() { sed -n 's/^#BSUB -o \(.*\)\/[^/]*$/\1/p' "$1" | head -1; }
is_queued() { bjobs -J "$1" 2>/dev/null | grep -qE "RUN|PEND|PSUSP|USUSP|SSUSP"; }
ACTIVE=" "    # jobs submitted in this run or still queued: later steps may wait for them
# submit an .lsf file unless its jobs are queued already; optional dependency on job $2
submit() {
  local f=$1 up=${2:-} n; n=$(jobname "$f"); mkdir -p "$(errdir "$f")"
  if is_queued "$n"; then stamp "queued already: $n"; ACTIVE="$ACTIVE$n "; return; fi
  if [ -n "$up" ] && [[ "$ACTIVE" == *" $up "* ]]; then bsub -w "ended($up)" < "$f" | sed "s/^/  $n: /"
  else bsub < "$f" | sed "s/^/  $n: /"; fi
  ACTIVE="$ACTIVE$n "
}
complete() { stamp "complete: $1"; }
ndone() { ls $@ 2>/dev/null | wc -l; }
SIM_LSF=$(ls Fig2/run_*.lsf Fig3/run_*.lsf Tab2/run_*.lsf SuppTab3/run_*.lsf SuppTab4/run_*.lsf SuppTab9/run_*.lsf SuppTab12/run_*.lsf)
APP_JOBS="prepare realdata medscc medscc_cap"
RD_DONE="Fig4/output/glgc_replication_*/DONE Fig4/output/ukb_pleiotropy_*/DONE Fig4/output/glgc_pleiotropy_*/DONE"
MD_DONE="Fig4/output/mediation_scc/*/DONE"

final_job() {
  # wait only for jobs that LSF still knows (finished jobs are forgotten after a while)
  local w="" n names
  names="$(for f in $SIM_LSF; do jobname "$f"; done) $APP_JOBS"
  for n in $names; do [ -n "$(bjobs -a -J "$n" 2>/dev/null)" ] && w="$w${w:+ && }ended($n)"; done
  mkdir -p logs
  bsub -J final_make -q "$QUEUE" -n 4 -W 12:00 -M 64 -R "span[hosts=1] rusage[mem=64]" \
       -o logs/final_make.out -e logs/final_make.err ${w:+-w "$w"} <<'JOB'
module load R/4.5.1 plink/1.90-beta bcftools/1.22.1 2>/dev/null
bash make_all.sh
find Data/raw -type f \( -name "*.gz" -o -name "*.bgz" -o -name "*.db" -o -name "*.gtf" \) -print0 2>/dev/null | sort -z | xargs -0 md5sum > logs/raw_md5.txt
JOB
}

case $mode in
status)
  Rscript -e 'here::i_am("run_all.sh"); source(here::here("SupportingCode","load_all.R")); print(sim_progress(), row.names = FALSE)'
  bjobs -w 2>/dev/null | awk 'NR == 1 || /final_make|prepare|realdata|medscc/'
  echo "application fits done: $(ndone $RD_DONE) of 13 (realdata), $(ndone $MD_DONE) of 15 (medscc)"
  ;;
final)
  final_job
  ;;
make)
  bash make_all.sh
  ;;
all)
  stamp "1. setup"
  if ! ls SupportingCode/src/reprokernels.* >/dev/null 2>&1 || [ ! -s logs/versions.txt ]; then
    Rscript setup.R > logs/setup.log 2>&1 || { stamp "setup failed, see logs/setup.log"; exit 1; }
  fi
  tail -n 12 logs/versions.txt

  stamp "2. simulation arrays"
  done_arrays=" $(Rscript -e 'here::i_am("run_all.sh"); source(here::here("SupportingCode","load_all.R")); p <- sim_progress(); cat(p$array[p$found >= p$expected])' 2>/dev/null | tail -1) "
  [ "$(ndone SuppTab4/output/floor_rep*.csv)" -ge 10 ] && done_arrays="$done_arrays SuppTab4 "
  for f in $SIM_LSF; do
    n=$(jobname "$f"); if [[ "$done_arrays" == *" $n "* ]]; then complete "$n"; else submit "$f"; fi
  done

  stamp "3. downloads"
  prepared=0
  [ "$(ndone Data/prepared/{glgc_replication,glgc_pleiotropy,ukb_pleiotropy}/DONE)" -eq 3 ] && prepared=1
  if [ $prepared = 1 ]; then complete "downloads (all data prepared)"
  elif [ "${DRYRUN:-0}" = "1" ]; then stamp "[dry run] would download the data"
  elif [ "${SKIP_DOWNLOAD:-0}" != "1" ]; then
    ok=1
    bash Data/download_data.sh > logs/download_data.log 2>&1 || ok=0
    grep -q "DONE" logs/download_data.log || { ok=0; stamp "incomplete: logs/download_data.log"; tail -3 logs/download_data.log; }
    [ $ok = 1 ] || { stamp "a download failed; rerun this script to resume"; exit 1; }
  fi

  stamp "4. data preparation and application fits"
  if [ "$(ndone Data/prepared/{glgc_replication,glgc_pleiotropy,ukb_pleiotropy}/DONE)" -eq 3 ]; then complete prepare; else submit Data/run_prepare.lsf; fi
  if [ "$(ndone $RD_DONE)" -ge 13 ]; then complete realdata; else submit Fig4/run_realdata.lsf prepare; fi
  if [ ! -s Data/mediation/mediation_scc_cis.tsv.gz ]; then stamp "lung SCC mediation data not prepared (not public, see Data/README.md): skipped"
  elif [ "$(ndone $MD_DONE)" -ge 15 ]; then complete medscc; complete medscc_cap
  else submit Fig4/run_mediation_scc.lsf; submit Fig4/run_mediation_scc_cap.lsf medscc; fi

  stamp "5. final job: figures and tables"
  is_queued final_make && stamp "queued already: final_make" || final_job
  stamp "submitted; when final_make has finished, read logs/make/status.txt"
  ;;
*) echo "usage: bash run_all.sh [all|status|final|make]"; exit 1 ;;
esac
