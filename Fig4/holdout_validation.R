here::i_am("Fig4/holdout_validation.R")
suppressPackageStartupMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
out <- here::here("Fig4", "output", "holdout_validation.tsv")
rows <- list()
for (analysis in names(RD_ANALYSES)) {
  holdout <- RD_ANALYSES[[analysis]]$holdout
  if (is.null(holdout)) next
  path <- here::here("Fig4", "output", paste0(analysis, "_main"), "discoveries_t1.tsv.gz")
  if (!file.exists(path)) next
  d <- fread(path)
  present <- intersect(holdout, names(d))
  rows[[analysis]] <- data.table(analysis, n_discoveries = nrow(d),
    holdout_coordinates = length(present), replicated_all = if (length(present)) sum(rowSums(d[, ..present] > 0) == length(present)) else NA_integer_)
}
fwrite(rbindlist(rows, fill = TRUE), out, sep = "\t")
