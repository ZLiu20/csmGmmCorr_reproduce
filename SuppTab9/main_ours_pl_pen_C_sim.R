# Simulation entry point for main_ours_pl_pen_C.
here::i_am("SuppTab9/main_ours_pl_pen_C_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("main_ours_pl_pen_C", as.integer(args[[1L]]))
