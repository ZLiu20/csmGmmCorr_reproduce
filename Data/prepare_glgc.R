# Prepare the GLGC analyses (SupportingCode/glgc_prepare.R; Supplementary Section 8.2).
#   Rscript Data/prepare_glgc.R replication     -> Data/prepared/glgc_replication/
#   Rscript Data/prepare_glgc.R pleiotropy      -> Data/prepared/glgc_pleiotropy/
# About 1-2 hours and 48 GB of memory per analysis (bsub < Data/run_prepare.lsf).
# Expected: replication 4,473,732 variants in all ten strata, 258,705 fitting variants;
#           pleiotropy 9,917,860 variants, 264,749 fitting variants.
here::i_am("Data/prepare_glgc.R")
source(here::here("SupportingCode", "load_all.R"))
an <- commandArgs(trailingOnly = TRUE)[1]
stopifnot(an %in% c("replication", "pleiotropy"))
prepare_glgc(manifest = here::here("Data", sprintf("glgc_%s_manifest.tsv", an)), raw_root = here::here("Data", "raw", "glgc"),
             out_dir = here::here("Data", "prepared", paste0("glgc_", an)), overwrite = TRUE)
