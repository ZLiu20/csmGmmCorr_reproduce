here::i_am("SuppTab6/make_supptab6_8.R")
source(here::here("SupportingCode", "load_all.R"))
d <- read_sim("Fig2", "main")
for (metric in c("FDP", "power", "R")) {
  number <- match(metric, c("FDP", "power", "R")) + 5L
  write_metric_table(d, metric, here::here("SuppTab6", sprintf("SuppTab%d.tex", number)),
                     sprintf("Design-level %s in the main simulation.", metric),
                     sprintf("tab:main-%s", tolower(metric)), c("target", "K", "t", "design", "method"))
}
