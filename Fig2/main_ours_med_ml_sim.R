# Fig2: main simulation (AR(1) correlation 0.2, logistic intercept -1); mediation;
# correlated csmGmm, maximum likelihood.
# Designs 1-12, K = 4 and 6, replicates 1-200; one task = 10 replicates of one
# (design, K, t) cell, 960 tasks. Output: Fig2/output/main/ours_med_ml/
#   Rscript Fig2/main_ours_med_ml_sim.R <task 1-960>        (or bsub < Fig2/run_main_ours_med_ml.lsf)
here::i_am("Fig2/main_ours_med_ml_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("main_ours_med_ml", aID)
