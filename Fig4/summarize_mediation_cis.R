# Summary of the two-mediator (blood, liver) mediation application.  Rscript Fig4/summarize_mediation_cis.R
here::i_am("Fig4/summarize_mediation_cis.R")
suppressMessages(library(data.table))
od <- here::here("Fig4", "output", "mediation_cis")
for (v in list.dirs(od, recursive = FALSE, full.names = FALSE)) {
  cat("\n==========", v, "==========\n")
  print(fread(file.path(od, v, "summary.tsv")))
  print(fread(file.path(od, v, "sigma.tsv"))[, .(stat1, stat2, fitted = round(fitted, 3), robust = round(robust, 3))])
  w <- fread(file.path(od, v, "support_weights.tsv")); print(w[, .(weight = round(sum(weight), 4)), by = support][order(-weight)])
  d <- fread(file.path(od, v, "discoveries.tsv"))
  cat("discoveries at t = 1 (first 40):\n")
  print(head(d[t == 1 & rej == TRUE][order(lfdr), .(gene_name, chr, alpha_blood = round(alpha_blood, 1), beta_blood = round(beta_blood, 1),
                                              alpha_liver = round(alpha_liver, 1), beta_liver = round(beta_liver, 1), lfdr = signif(lfdr, 2))], 40))
  cat("discoveries at t = 3:", paste(d[t == 3 & rej == TRUE]$gene_name, collapse = ", "), "\n")
}
