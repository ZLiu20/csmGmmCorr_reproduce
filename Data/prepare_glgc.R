here::i_am("Data/prepare_glgc.R")
source(here::here("SupportingCode", "load_all.R"))
analysis <- commandArgs(trailingOnly = TRUE)[1L]
stopifnot(analysis %in% c("replication", "pleiotropy"))
prepare_glgc(
  manifest = here::here("Data", sprintf("glgc_%s_manifest.tsv", analysis)),
  raw_root = here::here("Data", "raw", "glgc"),
  out_dir = here::here("Data", "prepared", paste0("glgc_", analysis)),
  overwrite = TRUE
)
