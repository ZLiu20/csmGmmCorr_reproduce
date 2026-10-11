# Simulation entry point for main_ours_rp_ml.
here::i_am("Fig2/main_ours_rp_ml_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("main_ours_rp_ml", as.integer(args[[1L]]))
