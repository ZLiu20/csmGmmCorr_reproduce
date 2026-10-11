here::i_am("SuppTab4/make_supptab4.R")
source(here::here("SupportingCode", "load_all.R"))
files <- list.files(here::here("SuppTab4", "output"), "[.]csv$", full.names = TRUE)
d <- data.table::rbindlist(lapply(files, data.table::fread), fill = TRUE)
write_metric_table(d, "FDP", here::here("SuppTab4", "SuppTab4.tex"),
                   "Sensitivity to the lower bound on support-size probabilities.",
                   "tab:weight-floor", c("fit", "target"))
