# Mediation application (main-text Section 6): one sentinel SNP, one fit, inference at t = 1, 2.
#   Rscript Fig4/mediation_fit.R <job 1-135>     (or bsub < Fig4/run_mediation.lsf)
# Jobs (SupportingCode/realdata.R, MED_JOBS): 15 SNPs x {estimated Sigma, Sigma = I, robust Sigma} x
# three tissue sets (blood/liver/adipose = primary; blood/adipose/aorta and blood/liver/aorta).
# Input: Data/mediation/mediation_scc_bysnp.tsv.gz. Output: Fig4/output/mediation/<tissues>/<snp>_<fit>/
here::i_am("Fig4/mediation_fit.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
j <- MED_JOBS[aID, ]
print(j)
run_mediation(j$snp, j$variant, j$tissues)
