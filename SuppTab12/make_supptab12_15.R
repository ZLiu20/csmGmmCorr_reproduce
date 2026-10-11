here::i_am("SuppTab12/make_supptab12_15.R")
source(here::here("SupportingCode", "load_all.R"))
d <- read_sim("SuppTab12", "main")
metrics <- c("FDP", "power", "R", "incon_score")
for (i in seq_along(metrics)) {
  write_metric_table(d, metrics[i], here::here("SuppTab12", sprintf("SuppTab%d.tex", i + 11L)),
    sprintf("Replication comparison: %s.", metrics[i]),
    sprintf("tab:replication-%s", tolower(metrics[i])))
}
