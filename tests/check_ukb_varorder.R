# Check that the *.varorder biomarker files of the Neale lab contain the same values as the files
# used in the article (the old names were removed from the public bucket in 2026).
#   Rscript tests/check_ukb_varorder.R <folder with the original files> <folder with the .varorder files>
library(data.table)
a <- commandArgs(TRUE); cols <- c("variant", "minor_AF", "low_confidence_variant", "n_complete_samples", "tstat")
for (f in c("30870_irnt", "30760_irnt", "30750_irnt")) {
  x <- fread(cmd = paste("gzip -cd", shQuote(file.path(a[1], paste0(f, ".gwas.imputed_v3.both_sexes.tsv.bgz")))), select = cols)
  y <- fread(cmd = paste("gzip -cd", shQuote(file.path(a[2], paste0(f, ".gwas.imputed_v3.both_sexes.varorder.tsv.bgz")))), select = cols)
  setkey(x, variant); setkey(y, variant)
  cat(sprintf("%s: %d vs %d rows, identical values: %s\n", f, nrow(x), nrow(y), isTRUE(all.equal(x, y, check.attributes = FALSE))))
}
