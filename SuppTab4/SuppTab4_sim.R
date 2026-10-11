here::i_am("SuppTab4/SuppTab4_sim.R")
source(here::here("SupportingCode", "load_all.R"))
replicate_id <- as.integer(commandArgs(trailingOnly = TRUE)[1L])
if (is.na(replicate_id)) replicate_id <- 1L
out <- floor_study_rep(replicate_id)
dir.create(here::here("SuppTab4", "output"), showWarnings = FALSE)
data.table::fwrite(out, here::here("SuppTab4", "output", sprintf("rep_%03d.csv", replicate_id)))
