#!/bin/bash
# Authors' check only: link the results of the original HPC runs into the folder layout of this
# repository, so that every plot_* and make_* script can be run without rerunning the
# simulations. Run from the repository root. Nothing is written into the original folders.
SIM=${SIM:-/rsrch8/home/biostatistics/zliu20/csmGmmCorr_SIM_20261005_v6001/results}
RD=${RD:-/rsrch8/home/biostatistics/zliu20/csmGmm_RealData/csmGmmCorr_REALDATA_v1}
PREP=${PREP:-/rsrch8/home/biostatistics/zliu20/csmGmm_RealData/formal_analysis/output/formal}
ln_dir() { [ -d "$1" ] || { echo "MISSING $1"; return; }; mkdir -p "$(dirname "$2")"; [ -e "$2" ] || ln -s "$1" "$2"; }
for sc in pl rp med; do
  for ty in pen ml; do
    ln_dir $SIM/main/ours_${sc}_$ty Fig2/output/main/ours_${sc}_$ty
    ln_dir $SIM/sym/ours_${sc}_$ty  Fig2/output/sym/ours_${sc}_$ty
    ln_dir $SIM/rho/ours_${sc}_$ty  Fig3/output/rho/ours_${sc}_$ty
  done
  ln_dir $SIM/rho/comp_${sc}_csm  Fig3/output/rho/comp_${sc}_csm
  ln_dir $SIM/main/comp_${sc}_csm Tab2/output/main/comp_${sc}_csm
  for st in quant exch split; do ln_dir $SIM/$st/ours_${sc}_pen Tab2/output/$st/ours_${sc}_pen; done
  for pg in L0A10 L10A10 L100A10 L30A0 L30A3 L30A30; do ln_dir $SIM/main/ours_${sc}_pen_$pg SuppTab3/output/main/ours_${sc}_pen_$pg; done
  for s0 in A C; do ln_dir $SIM/main/ours_${sc}_pen_$s0 SuppTab9/output/main/ours_${sc}_pen_$s0; done
done
ln_dir $SIM/ldu/ours_pl_pen Tab2/output/ld/ours_pl_pen     # LD study with all metrics (same data and fits as ld)
ln_dir $SIM/ldu/ours_rp_pen Tab2/output/ld/ours_rp_pen
ln_dir $SIM/ld/ours_med_pen Tab2/output/ld/ours_med_pen
ln_dir $SIM/rho/ours_med_pen_oracle Tab2/output/rho/ours_med_pen_oracle
ln_dir $SIM/main/comp_pl  Fig2/output/main/comp_pl
ln_dir $SIM/main/comp_med Fig2/output/main/comp_med
ln_dir $SIM/main/comp_rp  SuppTab12/output/main/comp_rp
# real data (results and prepared inputs)
for v in main indep robust main_amu15 indep_amu15; do
  ln_dir $RD/results/replication_$v Fig4/output/glgc_replication_$v
done
for v in main indep robust main_amu15; do
  ln_dir $RD/results/pleiotropy_ukb_$v Fig4/output/ukb_pleiotropy_$v
  ln_dir $RD/results/pleiotropy_$v Fig4/output/glgc_pleiotropy_$v
done
ln_dir $PREP/replication/prepared Data/prepared/glgc_replication
ln_dir $PREP/pleiotropy/prepared  Data/prepared/glgc_pleiotropy
ln_dir $RD/ukb/prepared            Data/prepared/ukb_pleiotropy
echo "links done"
