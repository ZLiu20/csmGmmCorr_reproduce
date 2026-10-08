# Tab2: quantitative outcomes; mediation;
# correlated csmGmm, penalized.
# Designs 1-7, K = 4 and 6, replicates 1-100; one task = 10 replicates of one
# (design, K, t) cell, 120 tasks. Output: Tab2/output/quant/ours_med_pen/
#   Rscript Tab2/quant_ours_med_pen_sim.R <task 1-120>        (or bsub < Tab2/run_quant_ours_med_pen.lsf)
here::i_am("Tab2/quant_ours_med_pen_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("quant_ours_med_pen", aID)
