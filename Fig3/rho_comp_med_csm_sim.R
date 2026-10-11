# Fig3: correlation study (rho = 0, 0.2, ..., 0.8); mediation;
# csmGmm (Sigma = I), plug-in Sigma and penalized csmGmm.
# Designs 1-5, K = 4 and 6, replicates 1-200; one task = 10 replicates of one
# (design, K, t) cell, 400 tasks. Output: Fig3/output/rho/comp_med_csm/
#   Rscript Fig3/rho_comp_med_csm_sim.R <task 1-400>        (or bsub < Fig3/run_rho_comp_med_csm.lsf)
here::i_am("Fig3/rho_comp_med_csm_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("rho_comp_med_csm", aID)
