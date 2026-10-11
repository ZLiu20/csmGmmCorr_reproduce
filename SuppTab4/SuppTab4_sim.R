# Supplementary Table 4: one replicate of the weight-floor study (Supplementary Section 2.6).
#   Rscript SuppTab4/SuppTab4_sim.R <replicate 1-10>      (or bsub < SuppTab4/run_SuppTab4.lsf)
here::i_am("SuppTab4/SuppTab4_sim.R")
source(here::here("SupportingCode", "load_all.R"))
aID <- as.integer(commandArgs(trailingOnly = TRUE)[1])
dir.create(here::here("SuppTab4", "output"), showWarnings = FALSE)
res <- floor_study_rep(aID)
print(res)
write.csv(res, here::here("SuppTab4", "output", sprintf("floor_rep%02d.csv", aID)), row.names = FALSE)
