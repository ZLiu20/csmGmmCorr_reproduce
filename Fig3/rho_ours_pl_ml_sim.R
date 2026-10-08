# Fig3: correlation study (rho = 0, 0.2, ..., 0.8); pleiotropy;
# correlated csmGmm, maximum likelihood.
# Designs 1-5, K = 4 and 6, replicates 1-200; one task = 10 replicates of one
# (design, K, t) cell, 200 tasks. Output: Fig3/output/rho/ours_pl_ml/
#   Rscript Fig3/rho_ours_pl_ml_sim.R <task 1-200>        (or bsub < Fig3/run_rho_ours_pl_ml.lsf)
here::i_am("Fig3/rho_ours_pl_ml_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("rho_ours_pl_ml", aID)
