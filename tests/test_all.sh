#!/usr/bin/env bash
set -euo pipefail
mode=${1:-local}
root="$(cd "$(dirname "$0")/.." && pwd)"
cd "$root"
mkdir -p tests/logs tests/errDir
case "$mode" in
  local)
    export VAL_N=${VAL_N:-1000}
    export VAL_TOL=${VAL_TOL:-1e-4}
    export VAL_MAPS=${VAL_MAPS:-5}
    export SIM_N=${SIM_N:-1000}
    export DIAG_MAX_ITER=${DIAG_MAX_ITER:-100}
    export INCON_MAX=${INCON_MAX:-1000}
    Rscript tests/check_syntax.R "$root"
    Rscript tests/compare_glgc_prep.R
    Rscript tests/test_prepare_mediation_ldl.R "$root"
    Rscript Validation/validate_csmgmm.R
    Rscript Validation/validate2_csmgmm.R
    Rscript Validation/diag_incon.R
    ;;
  report)
    Rscript -e 'here::i_am("setup.R"); source(here::here("SupportingCode","load_all.R")); print(sim_progress(), row.names=FALSE)'
    ;;
  submit)
    echo "Use 'DRYRUN=1 bash run_all.sh all' to inspect submissions, then run without DRYRUN."
    ;;
  *) echo "Usage: bash tests/test_all.sh {local|report|submit}" >&2; exit 2 ;;
esac
