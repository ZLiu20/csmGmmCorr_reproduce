here::i_am("SuppTab19/make_supptab19_23.R")
source(here::here("SupportingCode", "load_all.R"))
studies <- c("main", "split", "quant", "exch", "ld")
folders <- c("Fig2", rep("Tab2", 4L))
for (i in seq_along(studies)) {
  d <- read_sim(folders[i], studies[i])
  write_metric_table(d, "FDP", here::here("SuppTab19", sprintf("SuppTab%d.tex", i + 18L)),
    sprintf("Robustness design %s.", studies[i]), sprintf("tab:robust-%s", studies[i]))
}
