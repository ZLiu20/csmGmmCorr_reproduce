# Mediation application with three mediators per variant (liver expression of the three nearest
# genes; outcome LDL cholesterol; K = 6, t = 1, 2). Data: Data/prepare_mediation_ldl.R.
#   Rscript Fig4/mediation_ldl_fit.R <job 1-20>     (or bsub < Fig4/run_mediation_ldl.lsf)
# Jobs (MEDLDL_JOBS): data set "noshare" x {estimated Sigma, Sigma = I, robust Sigma} x a_mu in
# {0.5, 1.5} (jobs 1-6), data set "nearest" with estimated Sigma, a_mu in {0.5, 1.5} (7-8), and
# data set "noshare_dedup" (one variant per gene triple) x the three Sigma x a_mu in {0.5, 1.5} (9-14), and
# data set "noshare_cond" (alpha conditional on the lead eQTL variant) likewise (15-20).
here::i_am("Fig4/mediation_ldl_fit.R")
source(here::here("SupportingCode", "load_all.R"))
j <- MEDLDL_JOBS[as.integer(commandArgs(trailingOnly = TRUE)[1]), ]
print(j)
run_mediation_ldl(j$set, j$variant, a_mu = j$a_mu)
