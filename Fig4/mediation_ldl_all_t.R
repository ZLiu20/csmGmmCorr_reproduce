here::i_am("Fig4/mediation_ldl_all_t.R")
source(here::here("SupportingCode", "load_all.R"))

argument <- commandArgs(trailingOnly = TRUE)[1L]
out_root <- here::here("Fig4", "output", "mediation_ldl")
if (identical(argument, "summary")) {
  files <- file.path(medldl_fit_dirs(out_root), "all_t.tsv")
  files <- files[file.exists(files)]
  result <- data.table::rbindlist(lapply(files, data.table::fread), fill = TRUE)
  print(result, nrows = 500L)
} else {
  run_mediation_all_t(as.integer(argument), out_root = out_root)
}
