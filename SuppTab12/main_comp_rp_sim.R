# SuppTab12: main simulation (AR(1) correlation 0.2, logistic intercept -1); replication;
# repfdr and directional partial conjunction.
# Designs 1-12, K = 4 and 6, replicates 1-200; one task = 10 replicates of one
# (design, K, t) cell, 480 tasks. Output: SuppTab12/output/main/comp_rp/
#   Rscript SuppTab12/main_comp_rp_sim.R <task 1-480>        (or bsub < SuppTab12/run_main_comp_rp.lsf)
here::i_am("SuppTab12/main_comp_rp_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
run_sim_task("main_comp_rp", aID)
