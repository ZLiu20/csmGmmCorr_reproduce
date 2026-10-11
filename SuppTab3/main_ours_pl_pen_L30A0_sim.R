# SuppTab3: main simulation (AR(1) correlation 0.2, logistic intercept -1); pleiotropy;
# correlated csmGmm, penalized, lambda_mu = 30, lambda_w = 0.
# Designs 1-7, K = 4 and 6, replicates 1-50; one task = 10 replicates of one
# (design, K, t) cell, 30 tasks. Output: SuppTab3/output/main/ours_pl_pen_L30A0/
#   Rscript SuppTab3/main_ours_pl_pen_L30A0_sim.R <task 1-30>        (or bsub < SuppTab3/run_main_ours_pl_pen_L30A0.lsf)
here::i_am("SuppTab3/main_ours_pl_pen_L30A0_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("main_ours_pl_pen_L30A0", aID)
