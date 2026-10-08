# SuppTab3: main simulation (AR(1) correlation 0.2, logistic intercept -1); mediation;
# correlated csmGmm, penalized, lambda_mu = 10, lambda_w = 10.
# Designs 1-7, K = 4 and 6, replicates 1-50; one task = 10 replicates of one
# (design, K, t) cell, 60 tasks. Output: SuppTab3/output/main/ours_med_pen_L10A10/
#   Rscript SuppTab3/main_ours_med_pen_L10A10_sim.R <task 1-60>        (or bsub < SuppTab3/run_main_ours_med_pen_L10A10.lsf)
here::i_am("SuppTab3/main_ours_med_pen_L10A10_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("main_ours_med_pen_L10A10", aID)
