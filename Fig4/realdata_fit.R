# Applications (main-text Section 6): one fit and its inference at t = 1, 2.
#   Rscript Fig4/realdata_fit.R <job 1-13>     (or bsub < Fig4/run_realdata.lsf)
# Jobs (SupportingCode/realdata.R, RD_JOBS):
#    1-5   GLGC replication: estimated Sigma, Sigma = I, robust Sigma; a_mu = 1.5 with estimated Sigma and Sigma = I
#    6-9   UK Biobank pleiotropy: estimated Sigma, Sigma = I, robust Sigma; a_mu = 1.5
#   10-13  GLGC lipid pleiotropy (Supplementary Section S5.5): as 6-9
# Requires the prepared data in Data/prepared/ (Data/README.md). Output: Fig4/output/<analysis>_<fit>/
here::i_am("Fig4/realdata_fit.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
j <- RD_JOBS[aID, ]
print(j)
run_realdata(j$analysis, j$variant, a_mu = j$a_mu)
