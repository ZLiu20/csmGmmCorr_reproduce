# Mediation application (main-text Section 6): one sentinel SNP, one fit, inference at t = 1, 2.
#   Rscript Fig4/mediation_fit.R <job 1-855>     (or bsub < Fig4/run_mediation.lsf)
# Jobs (SupportingCode/realdata.R, MED_JOBS): 15 SNPs x {estimated Sigma, Sigma = I, robust Sigma} x
# tissue sets (three tissues: blood/liver/adipose, blood/adipose/aorta, blood/liver/aorta; one tissue:
# blood, liver, adipose, aorta, coronary artery) x lower bound a_mu of the nonzero magnitudes.
# Input: Data/mediation/mediation_scc_bysnp.tsv.gz. Output: Fig4/output/mediation/<tissues>/<snp>_<fit>/
here::i_am("Fig4/mediation_fit.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
j <- MED_JOBS[aID, ]
print(j)
run_mediation(j$snp, j$variant, j$tissues, a_mu = j$a_mu)
