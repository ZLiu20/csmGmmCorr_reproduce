# Simulation entry point for quant_ours_med_pen.
here::i_am("Tab2/quant_ours_med_pen_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("quant_ours_med_pen", as.integer(args[[1L]]))
