here::i_am("SuppTab3/make_supptab3.R")
suppressPackageStartupMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
penalty <- read_sim("SuppTab3", "main")
baseline <- three_designs(read_sim("Fig2", "main"))[method == "penalized"]
d <- rbind(penalty, baseline, fill = TRUE)
d <- d[rep <= 50L]
fwrite(d, here::here("SuppTab3", "SuppTab3.csv"))
write_metric_table(d, "FDP", here::here("SuppTab3", "SuppTab3.tex"),
                   "Sensitivity to the magnitude and weight penalty constants.",
                   "tab:penalty-choice", c("target", "K", "t", "design", "method"))
