# Simulation entry point for exch_ours_rp_pen.
here::i_am("Tab2/exch_ours_rp_pen_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("exch_ours_rp_pen", as.integer(args[[1L]]))
