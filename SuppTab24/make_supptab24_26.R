here::i_am("SuppTab24/make_supptab24_26.R")
source(here::here("SupportingCode", "load_all.R"))
d <- rd_collect(here::here("Fig4", "output", "mediation_ldl"))$summary
for (i in seq_along(c("R_raw", "R_adjusted", "incon_raw"))) {
  metric <- c("R_raw", "R_adjusted", "incon_raw")[i]
  if (!metric %in% names(d)) next
  write_metric_table(d, metric, here::here("SuppTab24", sprintf("SuppTab%d.tex", i + 23L)),
    sprintf("Mediation LDL analysis: %s.", metric), sprintf("tab:med-ldl-%s", tolower(metric)), c("variant", "t"))
}
