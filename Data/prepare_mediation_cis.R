# Mediation application (main-text Section 6; Figure of the two-mediator model): the exposure is the
# cis genetic variation of a gene and the two mediators are its expression in whole blood (M1) and
# liver (M2); the outcome is lung squamous cell carcinoma. One feature per gene (autosomal genes
# with all four statistics), K = 4, coordinates in the mediation order:
#   alpha_blood : exposure -> M1, gene-level cis-eQTL test in GTEx v8 whole blood: the
#                 permutation-calibrated p-value of the lead cis variant (pval_beta, which accounts for
#                 the selection of the lead variant), z = sign(slope) * qnorm(pval_beta / 2, lower = FALSE)
#   beta_blood  : M1 -> outcome, S-PrediXcan z of the gene's predicted whole-blood expression (MASHR)
#   alpha_liver, beta_liver : the same for liver.
#   Rscript Data/prepare_mediation_cis.R <folder with GTEx *.v8.egenes.txt.gz> <folder with scc__PM__*.csv>
here::i_am("Data/prepare_mediation_cis.R")
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); gtex <- a[1]; twas <- a[2]
tis <- c(blood = "Whole_Blood", liver = "Liver")
eg <- lapply(tis, function(t) {
  d <- fread(cmd = paste("gzip -cd", shQuote(file.path(gtex, sprintf("%s.v8.egenes.txt.gz", t)))))
  d <- d[is.finite(slope) & pval_beta > 0]
  data.table(gene_id = sub("\\..*$", "", d$gene_id), gene_name = d$gene_name, chr = sub("^chr", "", d$gene_chr),
             start_b38 = d$gene_start, z = sign(d$slope) * qnorm(d$pval_beta / 2, lower.tail = FALSE), lead = d$variant_id)
})
tw <- lapply(tis, function(t) {
  d <- fread(file.path(twas, sprintf("scc__PM__%s.csv", t)), select = c("gene", "zscore"))
  d <- d[is.finite(zscore)][, gene_id := sub("\\..*$", "", gene)]
  d[!(gene_id %in% gene_id[duplicated(gene_id)]), .(gene_id, zscore)]
})
D <- merge(eg$blood[, .(gene_id, gene_name, chr, start_b38, alpha_blood = z, lead_blood = lead)],
           tw$blood[, .(gene_id, beta_blood = zscore)], by = "gene_id")
D <- merge(D, eg$liver[, .(gene_id, alpha_liver = z, lead_liver = lead)], by = "gene_id")
D <- merge(D, tw$liver[, .(gene_id, beta_liver = zscore)], by = "gene_id")
D <- D[chr %in% as.character(1:22)]
setcolorder(D, c("gene_id", "gene_name", "chr", "start_b38", "alpha_blood", "beta_blood", "alpha_liver", "beta_liver",
                 "lead_blood", "lead_liver"))
fwrite(D, here::here("Data", "mediation", "mediation_cis_blood_liver.tsv.gz"), sep = "\t")
cat("genes:", nrow(D), "\n"); print(round(cor(D[, .(alpha_blood, beta_blood, alpha_liver, beta_liver)]), 3))
