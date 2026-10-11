# Simulation entry point for main_comp_pl.
here::i_am("Fig2/main_comp_pl_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("main_comp_pl", as.integer(args[[1L]]))
