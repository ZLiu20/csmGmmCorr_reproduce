# Authors' check: run the R preparation of the GLGC replication data (Data/prepare_glgc.R logic)
# into tests/glgc_replication_R/ and compare it with the prepared files used for the paper.
#   Rscript tests/compare_glgc_prep.R <raw GLGC folder> <prepared folder used for the paper>
here::i_am("tests/compare_glgc_prep.R")
source(here::here("SupportingCode", "load_all.R"))
suppressMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE); raw <- a[1]; ref <- a[2]
out <- here::here("tests", "glgc_replication_R")
t0 <- Sys.time()
prepare_glgc(here::here("Data", "glgc_replication_manifest.tsv"), raw, out, overwrite = TRUE)
cat("time:", format(Sys.time() - t0), "\n")
md5 <- function(f) sub(" .*", "", system(paste("gzip -cd", shQuote(f), "| md5sum"), intern = TRUE))
for (f in c("fit_sample.tsv.gz", "chr22.tsv.gz", "chr1.tsv.gz")) {
  a1 <- md5(file.path(out, f)); a2 <- md5(file.path(ref, f))
  cat(sprintf("%-18s R %s  reference %s  %s\n", f, a1, a2, if (a1 == a2) "IDENTICAL" else "DIFFERENT"))
}
m1 <- fread(file.path(out, "chromosome_manifest.tsv")); m2 <- fread(file.path(ref, "chromosome_manifest.tsv"))
cat("rows: R", sum(m1$n_rows), " reference", sum(m2$n_rows), "\n")
f1 <- fread(file.path(out, "fit_sample.tsv.gz"), select = "variant_id"); f2 <- fread(file.path(ref, "fit_sample.tsv.gz"), select = "variant_id")
cat("fitting variants: R", nrow(f1), " reference", nrow(f2), " common", length(intersect(f1$variant_id, f2$variant_id)), "\n")
cat("(the Mac copy used for the paper: fit_sample md5 3ec08daa6c6d75daf00e16c6b256b007, chr22 md5 f0ffb9d322f5b2e99c1da90488394905)\n")
