# Summary of the mediation application: discoveries (genes) per SNP, fit and tissue set at
# t = 1, 2, and the fitted correlations.   Rscript Fig4/summarize_mediation.R
here::i_am("Fig4/summarize_mediation.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
od <- here::here("Fig4", "output", "mediation")
S <- rbindlist(lapply(list.files(od, pattern = "^summary.tsv$", recursive = TRUE, full.names = TRUE), fread), fill = TRUE)
G <- rbindlist(lapply(list.files(od, pattern = "^sigma.tsv$", recursive = TRUE, full.names = TRUE), fread), fill = TRUE)
fwrite(S, file.path(od, "summary_mediation.tsv"), sep = "\t"); fwrite(G, file.path(od, "sigma_mediation.tsv"), sep = "\t")
cat("jobs done:", nrow(unique(S[, .(snp, variant, tissues)])), "of", nrow(MED_JOBS), "; not converged:", sum(!S$converged), "\n\n")
for (ts0 in names(MED_TISSUES)) {
  cat("==== tissues:", ts0, "====\n")
  w <- dcast(S[tissues == ts0], snp + snp_gene ~ variant + t, value.var = "R_raw")
  print(w)
  cat("totals:\n"); print(S[tissues == ts0, .(genes = sum(R_raw), snps_with_discoveries = sum(R_raw > 0), adj = sum(R_adj),
                                             extrap = sum(extrap_flags), null_weight = round(median(null_weight), 4)), by = .(variant, t)])
  cat("median fitted vs robust correlation (estimated Sigma):\n")
  print(G[tissues == ts0 & variant == "main", .(fitted = round(median(fitted), 3), robust = round(median(robust), 3)), by = .(stat1, stat2)][abs(robust) > 0.1])
  cat("\n")
}
cat("genes found with estimated Sigma, primary tissues, t = 1:\n")
print(S[tissues == names(MED_TISSUES)[1] & variant == "main" & t == 1 & R_raw > 0, .(snp, snp_gene, R_raw, genes = substr(genes, 1, 120))])
