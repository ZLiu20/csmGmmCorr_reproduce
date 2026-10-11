# Simulation entry point for main_ours_rp_pen_L30A0.
here::i_am("SuppTab3/main_ours_rp_pen_L30A0_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("main_ours_rp_pen_L30A0", as.integer(args[[1L]]))
