# Summary of the three-mediator (K = 6) LDL mediation application.  Rscript Fig4/summarize_mediation_ldl.R
here::i_am("Fig4/summarize_mediation_ldl.R")
source(here::here("SupportingCode", "load_all.R"))
suppressMessages(library(data.table))
od <- here::here("Fig4", "output", "mediation_ldl")
# genes of well-known LDL loci, used only to annotate the output
KNOWN <- c("SORT1", "CELSR2", "PSRC1", "PCSK9", "LDLR", "HMGCR", "APOE", "APOC1", "TOMM40", "APOB", "ABCG5", "ABCG8",
           "NPC1L1", "TRIB1", "HNF1A", "FADS1", "FADS2", "ANGPTL3", "LIPC", "MAFB", "TM6SF2", "SUGP1", "CETP", "LPA")
for (v in list.dirs(od, recursive = FALSE, full.names = FALSE)) {
  if (!file.exists(file.path(od, v, "summary.tsv"))) { cat("\n", v, ": no summary\n"); next }
  cat("\n==========", v, "==========\n")
  print(fread(file.path(od, v, "summary.tsv")))
  print(fread(file.path(od, v, "sigma.tsv"))[, .(stat1, stat2, fitted = round(fitted, 3), robust = round(robust, 3))])
  print(fread(file.path(od, v, "iqr_sd.tsv")))
  w <- fread(file.path(od, v, "support_weights.tsv")); print(w[, .(weight = round(sum(weight), 4)), by = support][order(-weight)])
  d <- fread(file.path(od, v, "discoveries.tsv.gz"))
  for (t0 in sort(unique(d$t))) {
    x <- d[t == t0 & rej == TRUE]
    if (!nrow(x)) next
    lead <- x[clump_leads(as.character(x$chr), x$pos_b38, x$lfdr, rep(TRUE, nrow(x)))]
    cat(sprintf("t = %d: %d variants, %d gene triples, %d loci; loci (lead variant, first 40):\n", t0, nrow(x), uniqueN(x[, .(gene1, gene2, gene3)]), nrow(lead)))
    print(head(lead[order(lfdr), .(chr, pos_b38, gene1_name, gene2_name, gene3_name, a1 = round(alpha1, 1), b1 = round(beta1, 1),
                                   a2 = round(alpha2, 1), b2 = round(beta2, 1), a3 = round(alpha3, 1), b3 = round(beta3, 1),
                                   lfdr = signif(lfdr, 2))], 40))
    k <- intersect(KNOWN, unlist(x[, .(gene1_name, gene2_name, gene3_name)]))
    cat("known LDL genes among the mediators:", paste(sort(k), collapse = ", "), "\n")
  }
}
