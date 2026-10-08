# Mediation application (main-text Section 6): statistics for the 15 lung-cancer sentinel SNPs of
# Sun, McCaw and Lin (2025, Section 6), one analysis per SNP. Feature = gene; for tissue g
# (whole blood, liver, subcutaneous adipose; GTEx v8):
#   alpha_g: SNP -> expression of the gene in tissue g (eQTL test statistic of the SNP for every
#            gene, from GTEx individual-level data as in Sun et al., Fig4/calculate_eqtl_stats);
#   beta_g : expression of the gene in tissue g -> lung squamous cell carcinoma (S-PrediXcan z,
#            PredictDB MASHR models).
# All genes with at least one statistic (missing values allowed); an analysis uses the genes with
# all statistics of its three tissues (SupportingCode/realdata.R, MED_TISSUES).
# The prepared file Data/mediation/mediation_scc_bysnp.tsv.gz is included in the repository
# (GTEx individual-level data cannot be shared); this script documents how it was built from
#   <raw>/rs<id>_<gene>_<chr>_<pos>_<tissue>.txt   (columns Gene, SNP, MAF, testStat, pval)
#   <raw>/scc__PM__<tissue>.csv                      (S-PrediXcan output; columns gene_name, zscore)
#   Rscript Data/prepare_mediation.R <raw folder>
here::i_am("Data/prepare_mediation.R")
suppressMessages(library(data.table))
raw <- commandArgs(trailingOnly = TRUE)[1]
tis <- c(blood = "Whole_Blood", liver = "Liver", adipose = "Adipose_Subcutaneous", aorta = "Artery_Aorta",
         coronary = "Artery_Coronary")
uniq <- function(d, key) d[!(get(key) %in% d[[key]][duplicated(d[[key]])])]   # drop names that occur twice
tw <- lapply(tis, function(t) { d <- fread(file.path(raw, sprintf("scc__PM__%s.csv", t)), select = c("gene_name", "zscore"))
  uniq(d[is.finite(zscore)], "gene_name") })
files <- list.files(raw, pattern = "^rs.*_Whole_Blood\\.txt$")
snps <- sub("_.*", "", files)
out <- rbindlist(lapply(seq_along(snps), function(i) {
  parts <- strsplit(files[i], "_")[[1]]
  d <- NULL
  for (k in names(tis)) {
    f <- list.files(raw, pattern = sprintf("^%s_.*_%s\\.txt$", snps[i], tis[[k]]), full.names = TRUE)
    e <- fread(f, select = c("Gene", "testStat")); e <- uniq(e[is.finite(testStat)], "Gene")
    x <- merge(e, tw[[k]], by.x = "Gene", by.y = "gene_name", all = TRUE)
    setnames(x, c("Gene", "testStat", "zscore"), c("gene", paste0("alpha_", k), paste0("beta_", k)))
    d <- if (is.null(d)) x else merge(d, x, by = "gene", all = TRUE)
  }
  cbind(data.table(snp = snps[i], snp_gene = parts[2], chr = parts[3], pos = as.integer(parts[4])), d)
}))
cols <- as.vector(outer(c("alpha_", "beta_"), names(tis), paste0))   # alpha_blood, beta_blood, alpha_liver, ...
setcolorder(out, c("snp", "snp_gene", "chr", "pos", "gene", cols))
fwrite(out, here::here("Data", "mediation", "mediation_scc_bysnp.tsv.gz"), sep = "\t")
print(out[, .N, by = snp])
