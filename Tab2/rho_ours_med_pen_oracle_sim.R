# Tab2: correlation study (rho = 0, 0.2, ..., 0.8); mediation;
# correlated csmGmm, penalized, Sigma fixed at the correlation of the null statistics.
# Designs 4-5, K = 4 and 6, replicates 1-50; one task = 10 replicates of one
# (design, K, t) cell, 40 tasks. Output: Tab2/output/rho/ours_med_pen_oracle/
#   Rscript Tab2/rho_ours_med_pen_oracle_sim.R <task 1-40>        (or bsub < Tab2/run_rho_ours_med_pen_oracle.lsf)
here::i_am("Tab2/rho_ours_med_pen_oracle_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("rho_ours_med_pen_oracle", aID)
