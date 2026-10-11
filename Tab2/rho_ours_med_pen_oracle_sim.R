# Simulation entry point for rho_ours_med_pen_oracle.
here::i_am("Tab2/rho_ours_med_pen_oracle_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("rho_ours_med_pen_oracle", as.integer(args[[1L]]))
