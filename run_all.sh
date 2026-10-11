#!/usr/bin/env bash
set -u
cd "$(dirname "$0")"
mkdir -p logs

QUEUE=${QUEUE:-medium}
DRYRUN=${DRYRUN:-0}
mode=${1:-all}
ACTIVE_IDS=()

stamp() { printf '[%s] %s\n' "$(date '+%F %T')" "$*"; }

job_name() {
  sed -n 's/^#BSUB[[:space:]]*-J[[:space:]]*//p' "$1" | head -1 |
    sed 's/\[[^]]*\].*$//'
}

active_job() {
  local name=$1
  bjobs -J "$name" -noheader -o jobid 2>/dev/null | head -1
}

submit_lsf() {
  local file=$1 dependency=${2:-}
  local name existing output id
  name=$(job_name "$file")
  existing=$(active_job "$name")
  if [[ -n "$existing" ]]; then
    stamp "already queued: $name ($existing)" >&2
    printf '%s\n' "$existing"
    return
  fi
  if [[ "$DRYRUN" == 1 ]]; then
    stamp "[dry run] $file${dependency:+ after $dependency}" >&2
    printf 'DRY_%s\n' "$name"
    return
  fi
  if [[ -n "$dependency" ]]; then
    output=$(bsub -w "ended($dependency)" < "$file") || return 1
  else
    output=$(bsub < "$file") || return 1
  fi
  id=$(printf '%s\n' "$output" | sed -n 's/.*Job <\([0-9][0-9]*\)>.*/\1/p' | head -1)
  [[ -n "$id" ]] || { printf '%s\n' "$output" >&2; return 1; }
  stamp "submitted $name as $id" >&2
  printf '%s\n' "$id"
}

complete() { stamp "complete: $*"; }

count_matching() {
  find "$1" -type f -name "$2" 2>/dev/null | wc -l | tr -d ' '
}

count_realdata_fits() {
  Rscript -e '
    here::i_am("run_all.sh")
    source(here::here("SupportingCode", "load_all.R"))
    paths <- rd_fit_dirs(here::here("Fig4", "output"))
    cat(sum(file.exists(file.path(paths, "summary.tsv"))))
  ' 2>/dev/null || printf '0'
}

status_report() {
  module load R/4.5.1 2>/dev/null || true
  Rscript -e 'here::i_am("run_all.sh"); source(here::here("SupportingCode", "load_all.R")); print(sim_progress(), row.names = FALSE)'
  bjobs -w 2>/dev/null || true
  stamp "prepared GLGC analyses: $(count_matching Data/prepared PREPARATION_COMPLETE)/3"
  stamp "real-data fits: $(count_realdata_fits)/13"
  stamp "mediation fits: $(count_matching Fig4/output/mediation_ldl summary.tsv)/23"
}

final_job() {
  local condition="" id output
  for id in "${ACTIVE_IDS[@]}"; do
    [[ "$id" == DRY_* ]] && continue
    condition+="${condition:+ && }ended($id)"
  done
  if [[ "$DRYRUN" == 1 ]]; then
    stamp "[dry run] final make${condition:+ after $condition}"
    return
  fi
  if [[ -n "$(active_job final_make)" ]]; then
    stamp "already queued: final_make"
    return
  fi
  local args=(-J final_make -q "$QUEUE" -W 12:00 -n 4 -M 64
              -R 'span[hosts=1] rusage[mem=64]'
              -o logs/final_make.out -e logs/final_make.err)
  [[ -n "$condition" ]] && args+=(-w "$condition")
  output=$(bsub "${args[@]}" bash make_all.sh) || return 1
  stamp "$output"
}

case "$mode" in
  status)
    status_report
    ;;
  make)
    bash make_all.sh
    ;;
  final)
    final_job
    ;;
  all)
    stamp "1. setup"
    if [[ ! -s logs/versions.txt ]] ||
       ! ls SupportingCode/src/reprokernels.* >/dev/null 2>&1; then
      Rscript setup.R >logs/setup.log 2>&1 || {
        tail -n 40 logs/setup.log
        exit 1
      }
    else
      complete setup
    fi

    stamp "2. simulation arrays"
    completed_arrays=$(Rscript -e '
      here::i_am("run_all.sh")
      source(here::here("SupportingCode", "load_all.R"))
      p <- sim_progress()
      cat(p$array[p$found >= p$expected], sep = " ")
    ' 2>/dev/null || true)
    if [[ "$(count_matching SuppTab4/output 'floor_rep*.csv')" -ge 10 ]]; then
      completed_arrays="$completed_arrays SuppTab4"
    fi
    for file in Fig2/run_*.lsf Fig3/run_*.lsf Tab2/run_*.lsf \
                SuppTab3/run_*.lsf SuppTab4/run_*.lsf \
                SuppTab9/run_*.lsf SuppTab12/run_*.lsf; do
      [[ -f "$file" ]] || continue
      name=$(job_name "$file")
      if [[ " $completed_arrays " == *" $name "* ]]; then
        complete "$name"
      else
        id=$(submit_lsf "$file") || exit 1
        ACTIVE_IDS+=("$id")
      fi
    done

    stamp "3. downloads"
    if [[ "${SKIP_DOWNLOAD:-0}" == 1 ]]; then
      stamp "downloads skipped by SKIP_DOWNLOAD"
    elif [[ "$DRYRUN" == 1 ]]; then
      stamp "[dry run] would download source data"
    else
      bash Data/download_data.sh >logs/download_data.log 2>&1 & p1=$!
      bash Data/download_mediation_ldl.sh >logs/download_mediation.log 2>&1 & p2=$!
      bash Data/download_1000g.sh >logs/download_1000g.log 2>&1 & p3=$!
      ok=1
      for process in "$p1" "$p2" "$p3"; do wait "$process" || ok=0; done
      [[ "$ok" == 1 ]] || {
        stamp "a download failed; inspect logs/download_*.log"
        exit 1
      }
    fi

    stamp "4. preparation and application fits"
    if [[ "$(count_matching Data/prepared PREPARATION_COMPLETE)" -ge 3 ]]; then
      complete prepare
      prepare_id=""
    else
      prepare_id=$(submit_lsf Data/run_prepare.lsf) || exit 1
      ACTIVE_IDS+=("$prepare_id")
    fi

    if [[ "$(count_realdata_fits)" -ge 13 ]]; then
      complete realdata
    else
      realdata_id=$(submit_lsf Fig4/run_realdata.lsf "$prepare_id") || exit 1
      ACTIVE_IDS+=("$realdata_id")
    fi

    if [[ -s Data/mediation/mediation_ldl_liver_noshare.tsv.gz ]]; then
      complete prep_mediation_ldl
      med_prepare_id=""
    else
      med_prepare_id=$(submit_lsf Data/run_prepare_mediation_ldl.lsf) || exit 1
      ACTIVE_IDS+=("$med_prepare_id")
    fi

    if [[ -s Data/mediation/mediation_ldl_liver_noshare_cond.tsv.gz ]]; then
      complete prep_mediation_ldl_cond
      med_cond_id="$med_prepare_id"
    else
      med_cond_id=$(submit_lsf Data/run_prepare_mediation_ldl_cond.lsf "$med_prepare_id") || exit 1
      ACTIVE_IDS+=("$med_cond_id")
    fi

    med_done=$(count_matching Fig4/output/mediation_ldl summary.tsv)
    if [[ "$med_done" -ge 23 ]]; then
      complete mediation_ldl
      med_fit_id="$med_cond_id"
    else
      med_fit_id=$(submit_lsf Fig4/run_mediation_ldl.lsf "$med_cond_id") || exit 1
      ACTIVE_IDS+=("$med_fit_id")
    fi

    all_t_done=$(count_matching Fig4/output/mediation_ldl all_t.tsv)
    if [[ "$all_t_done" -ge 23 ]]; then
      complete mediation_all_t
    else
      id=$(submit_lsf Fig4/run_mediation_ldl_all_t.lsf "$med_fit_id") || exit 1
      ACTIVE_IDS+=("$id")
    fi

    if [[ -s Fig4/output/mediation_ldl/method_comparison.tsv ]]; then
      complete mediation_compare
    else
      id=$(submit_lsf Fig4/run_mediation_ldl_compare.lsf "$med_fit_id") || exit 1
      ACTIVE_IDS+=("$id")
    fi

    stamp "5. final figures and tables"
    final_job
    ;;
  *)
    echo "usage: bash run_all.sh [all|status|final|make]" >&2
    exit 1
    ;;
esac
