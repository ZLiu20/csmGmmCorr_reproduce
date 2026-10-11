here::i_am("Fig4/mediation_ldl_compare.R")
suppressPackageStartupMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
sets <- c("noshare", "nearest", "noshare_dedup", "noshare_cond")
rows <- lapply(sets, function(set) {
  path <- here::here("Fig4", "output", "mediation_ldl", paste0(set, "_main"), "summary.tsv")
  if (file.exists(path)) fread(path) else NULL
})
d <- rbindlist(rows, fill = TRUE)
dir.create(here::here("Fig4", "output", "mediation_ldl"), recursive = TRUE, showWarnings = FALSE)
fwrite(d, here::here("Fig4", "output", "mediation_ldl", "method_comparison.tsv"), sep = "\t")
