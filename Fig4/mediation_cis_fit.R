# Mediation application with two mediators per feature (whole-blood and liver expression of a
# gene; exposure = the gene's cis genetic variation; outcome = lung squamous cell carcinoma).
#   Rscript Fig4/mediation_cis_fit.R <job 1-6>     (or bsub < Fig4/run_mediation_cis.lsf)
# Jobs (MEDCIS_JOBS): {estimated Sigma, Sigma = I, robust Sigma} x a_mu in {0.5, 1.5}.
here::i_am("Fig4/mediation_cis_fit.R")
source(here::here("SupportingCode", "load_all.R"))
j <- MEDCIS_JOBS[as.integer(commandArgs(trailingOnly = TRUE)[1]), ]
print(j)
run_mediation_cis(j$variant, a_mu = j$a_mu)
