# Prepare the UK Biobank pleiotropy analysis (Supplementary Section 8.1) in the format of the
# GLGC analyses: Data/prepared/ukb_pleiotropy/chr<k>.tsv.gz (variant_id, chr, pos_b37, BMI, SBP,
# TG, HDL, HbA1c), chromosome_manifest.tsv and fit_sample.tsv.gz (one variant per occupied
# 10-kb window, chosen at random with a fixed seed).
# Variants: autosomes, not low-confidence in any trait, minor allele frequency >= 0.01;
# statistic: tstat. Expected: 9,370,133 variants, 264,709 fitting variants.
#   Rscript Data/prepare_ukb.R      (about 30-60 min, 32 GB)
here::i_am("Data/prepare_ukb.R")
suppressMessages(library(data.table))
raw <- here::here("Data", "raw", "ukb"); out <- here::here("Data", "prepared", "ukb_pleiotropy")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
tr <- c(BMI = "21001_irnt", SBP = "4080_irnt", TG = "30870_irnt", HDL = "30760_irnt", HbA1c = "30750_irnt")
D <- NULL
for (k in names(tr)) {
  f <- file.path(raw, paste0(tr[[k]], ".gwas.imputed_v3.both_sexes.tsv.bgz"))
  x <- fread(cmd = paste("gzip -cd", shQuote(f)), select = c("variant", "minor_AF", "low_confidence_variant", "n_complete_samples", "tstat"))
  cat(sprintf("%-6s %s: %d variants, n = %d\n", k, basename(f), nrow(x), max(x$n_complete_samples, na.rm = TRUE)))
  x <- x[low_confidence_variant == FALSE & minor_AF >= 0.01 & is.finite(tstat), .(variant, maf = minor_AF, z = tstat)]
  setnames(x, "z", k)
  D <- if (is.null(D)) x else merge(D, x[, c("variant", k), with = FALSE], by = "variant")
  rm(x); gc(verbose = FALSE)
}
D[, c("chr", "pos_b37") := tstrsplit(variant, ":", keep = 1:2)]
D <- D[chr %in% as.character(1:22)]; D[, pos_b37 := as.integer(pos_b37)]
setnames(D, "variant", "variant_id"); setorder(D, chr, pos_b37)
cat("variants in all five traits:", nrow(D), "\n")
cols <- c("variant_id", "chr", "pos_b37", names(tr))
set.seed(20261007)
D[, win := paste(chr, pos_b37 %/% 10000L)]
fs <- D[D[, .I[sample(.N, 1L)], by = win]$V1, cols, with = FALSE]
fwrite(fs, file.path(out, "fit_sample.tsv.gz"), sep = "\t")
man <- rbindlist(lapply(as.character(1:22), function(ch) {
  p <- file.path(out, sprintf("chr%s.tsv.gz", ch)); y <- D[chr == ch, cols, with = FALSE]
  fwrite(y, p, sep = "\t"); data.table(chr = ch, path = basename(p), n_rows = nrow(y))
}))
fwrite(man, file.path(out, "chromosome_manifest.tsv"), sep = "\t")
Z <- as.matrix(D[, names(tr), with = FALSE])
cat("fit sample:", nrow(fs), "; test variants:", nrow(D) - nrow(fs), "\n")
cat("median-based lambda:", round(apply(Z, 2, function(z) median(z^2) / qchisq(0.5, 1)), 3), "\n")
cat("Pearson correlation (|z| < 2 in all five):\n"); print(round(cor(Z[rowSums(abs(Z) < 2) == ncol(Z), ]), 3))
writeLines("ok", file.path(out, "PREPARATION_COMPLETE"))
