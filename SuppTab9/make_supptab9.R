here::i_am("SuppTab9/make_supptab9.R")
source(here::here("SupportingCode", "load_all.R"))
d <- read_sim("SuppTab9", "main")
write_metric_table(d, "FDP", here::here("SuppTab9", "SuppTab9.tex"),
                   "Sensitivity to alternative initialization choices.",
                   "tab:initialization", c("target", "K", "t", "design", "method"))
