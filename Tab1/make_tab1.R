here::i_am("Tab1/make_tab1.R")
suppressPackageStartupMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
D <- design_means(read_sim("Fig2", "main"))
O <- overall_means(D)
write_metric_table(O, "FDP", here::here("Tab1", "Tab1.tex"),
                   "Simulation results averaged over the twelve designs of each scientific target.",
                   "tab:simulation-overview")
fwrite(O, here::here("Tab1", "Tab1_overall_means.csv"))
