# Simulation entry point for main_comp_med.
here::i_am("Fig2/main_comp_med_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("main_comp_med", as.integer(args[[1L]]))
