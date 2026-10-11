# Simulation entry point for main_ours_med_pen_A.
here::i_am("SuppTab9/main_ours_med_pen_A_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("main_ours_med_pen_A", as.integer(args[[1L]]))
