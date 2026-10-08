# Tab2: main simulation (AR(1) correlation 0.2, logistic intercept -1); pleiotropy;
# csmGmm (Sigma = I), plug-in Sigma and penalized csmGmm.
# Designs 1-9, K = 4 and 6, replicates 1-200; one task = 10 replicates of one
# (design, K, t) cell, 360 tasks. Output: Tab2/output/main/comp_pl_csm/
#   Rscript Tab2/main_comp_pl_csm_sim.R <task 1-360>        (or bsub < Tab2/run_main_comp_pl_csm.lsf)
here::i_am("Tab2/main_comp_pl_csm_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("main_comp_pl_csm", aID)
