# Simulation entry point for main_comp_med_csm.
here::i_am("Tab2/main_comp_med_csm_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("main_comp_med_csm", as.integer(args[[1L]]))
