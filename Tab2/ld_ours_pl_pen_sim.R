# Tab2: linkage disequilibrium between neighbouring variants; pleiotropy;
# correlated csmGmm, penalized.
# Designs 1-7, K = 4 and 6, replicates 1-100; one task = 10 replicates of one
# (design, K, t) cell, 60 tasks. Output: Tab2/output/ld/ours_pl_pen/
#   Rscript Tab2/ld_ours_pl_pen_sim.R <task 1-60>        (or bsub < Tab2/run_ld_ours_pl_pen.lsf)
here::i_am("Tab2/ld_ours_pl_pen_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("ld_ours_pl_pen", aID)
