# Fig2: zero-intercept logistic designs (logistic intercept 0); replication;
# correlated csmGmm, penalized.
# Designs 1-12, K = 4 and 6, replicates 1-200; one task = 10 replicates of one
# (design, K, t) cell, 480 tasks. Output: Fig2/output/sym/ours_rp_pen/
#   Rscript Fig2/sym_ours_rp_pen_sim.R <task 1-480>        (or bsub < Fig2/run_sym_ours_rp_pen.lsf)
here::i_am("Fig2/sym_ours_rp_pen_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("sym_ours_rp_pen", aID)
