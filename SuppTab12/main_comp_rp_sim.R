# Simulation entry point for main_comp_rp.
here::i_am("SuppTab12/main_comp_rp_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("main_comp_rp", as.integer(args[[1L]]))
