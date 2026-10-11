# SuppTab9: main simulation (AR(1) correlation 0.2, logistic intercept -1); mediation;
# correlated csmGmm, penalized, magnitude structure C (2/2).
# Designs 1-7, K = 4 and 6, replicates 1-50; one task = 10 replicates of one
# (design, K, t) cell, 60 tasks. Output: SuppTab9/output/main/ours_med_pen_C/
#   Rscript SuppTab9/main_ours_med_pen_C_sim.R <task 1-60>        (or bsub < SuppTab9/run_main_ours_med_pen_C.lsf)
here::i_am("SuppTab9/main_ours_med_pen_C_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("main_ours_med_pen_C", aID)
