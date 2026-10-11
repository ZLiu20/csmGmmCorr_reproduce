here::i_am("SuppTab16/make_supptab16_18.R")
source(here::here("SupportingCode", "load_all.R"))
d <- read_sim("Fig3", "rho")
for (i in seq_along(c("FDP", "power", "R"))) {
  metric <- c("FDP", "power", "R")[i]
  write_metric_table(d, metric, here::here("SuppTab16", sprintf("SuppTab%d.tex", i + 15L)),
    sprintf("Correlation sensitivity: %s.", metric), sprintf("tab:rho-%s", tolower(metric)))
}
