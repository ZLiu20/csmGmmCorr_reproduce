# Fig2: main simulation (AR(1) correlation 0.2, logistic intercept -1); mediation;
# correlated csmGmm, penalized.
# Designs 1-12, K = 4 and 6, replicates 1-500; one task = 10 replicates of one
# (design, K, t) cell, 2400 tasks. Output: Fig2/output/main/ours_med_pen/
#   Rscript Fig2/main_ours_med_pen_sim.R <task 1-2400>        (or bsub < Fig2/run_main_ours_med_pen.lsf)
here::i_am("Fig2/main_ours_med_pen_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("main_ours_med_pen", aID)
