here::i_am("Fig4/diagnose_mediation_ldl.R")
suppressPackageStartupMessages(library(data.table))
sets <- c("noshare", "nearest", "noshare_dedup", "noshare_cond")
rows <- lapply(sets, function(set) {
  path <- here::here("Data", "mediation", sprintf("mediation_ldl_liver_%s.tsv.gz", set))
  if (!file.exists(path)) return(NULL)
  d <- fread(path)
  zcols <- grep("^(alpha|beta)[0-9]+$", names(d), value = TRUE)
  data.table(set, n = nrow(d), n_gene_triples = uniqueN(do.call(paste, c(d[, grep("^gene[0-9]+$", names(d)), with=FALSE], sep="|"))),
             max_abs_z = max(abs(as.matrix(d[, ..zcols])), na.rm = TRUE))
})
fwrite(rbindlist(rows, fill = TRUE), here::here("Fig4", "mediation_ldl_diagnostics.tsv"), sep = "\t")
