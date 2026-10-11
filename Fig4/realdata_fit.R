here::i_am("Fig4/realdata_fit.R")
source(here::here("SupportingCode", "load_all.R"))
index <- as.integer(commandArgs(trailingOnly = TRUE)[1L])
if (is.na(index) || index < 1L || index > nrow(RD_JOBS)) stop("Supply a valid RD_JOBS row index.")
job <- RD_JOBS[index, ]
run_realdata(job$analysis, job$variant, job$a_mu)
