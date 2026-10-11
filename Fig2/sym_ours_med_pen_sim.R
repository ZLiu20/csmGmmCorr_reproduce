# Simulation entry point for sym_ours_med_pen.
here::i_am("Fig2/sym_ours_med_pen_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("sym_ours_med_pen", as.integer(args[[1L]]))
