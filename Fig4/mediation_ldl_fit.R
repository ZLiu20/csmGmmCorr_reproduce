here::i_am("Fig4/mediation_ldl_fit.R")
source(here::here("SupportingCode", "load_all.R"))
args <- commandArgs(trailingOnly = TRUE)
set <- if (length(args)) args[1L] else "noshare"
variant <- if (length(args) >= 2L) args[2L] else "main"
a_mu <- if (length(args) >= 3L) as.numeric(args[3L]) else 0.5
run_mediation_ldl(set, variant, a_mu)
