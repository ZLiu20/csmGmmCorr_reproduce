here::i_am("Data/prepare_ukb.R")
suppressPackageStartupMessages(library(data.table))
raw <- here::here("Data", "raw", "ukb")
out <- here::here("Data", "prepared", "ukb_pleiotropy")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
traits <- c(BMI = "21001_irnt", SBP = "4080_irnt", TG = "30870_irnt",
            HDL = "30760_irnt", HbA1c = "30750_irnt")
merged <- NULL
for (label in names(traits)) {
  path <- file.path(raw, paste0(traits[[label]], ".gwas.imputed_v3.both_sexes.tsv.bgz"))
  x <- fread(cmd = paste("gzip -cd", shQuote(path)),
             select = c("variant", "minor_AF", "low_confidence_variant", "n_complete_samples", "tstat"))
  message(sprintf("%-6s %s: %d variants", label, basename(path), nrow(x)))
  x <- x[low_confidence_variant == FALSE & minor_AF >= 0.01 & is.finite(tstat), .(variant, z = tstat)]
  setnames(x, "z", label)
  merged <- if (is.null(merged)) x else merge(merged, x, by = "variant")
}
loc <- tstrsplit(merged$variant, ":", fixed = TRUE, keep = 1:2)
merged[, `:=`(chr = as.character(loc[[1L]]), pos_b37 = as.integer(loc[[2L]]), variant_id = variant)]
setorder(merged, chr, pos_b37)
cols <- c("variant_id", "chr", "pos_b37", names(traits))
set.seed(20261007)
merged[, window := paste(chr, pos_b37 %/% 10000L)]
fit <- merged[, .SD[sample.int(.N, 1L)], by = window, .SDcols = cols]
fwrite(fit[, ..cols], file.path(out, "fit_sample.tsv.gz"), sep = "\t")
manifest <- rbindlist(lapply(as.character(1:22), function(chromosome) {
  y <- merged[chr == chromosome, ..cols]
  path <- file.path(out, sprintf("chr%s.tsv.gz", chromosome))
  fwrite(y, path, sep = "\t")
  data.table(chr = chromosome, path = basename(path), n_rows = nrow(y))
}))
fwrite(manifest, file.path(out, "chromosome_manifest.tsv"), sep = "\t")
writeLines("ok", file.path(out, "PREPARATION_COMPLETE"))
