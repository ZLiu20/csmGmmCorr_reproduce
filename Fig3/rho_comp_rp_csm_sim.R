# Simulation entry point for rho_comp_rp_csm.
here::i_am("Fig3/rho_comp_rp_csm_sim.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
if (!length(args)) stop("An LSF task index is required")
run_sim_task("rho_comp_rp_csm", as.integer(args[[1L]]))
