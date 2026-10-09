# One-time setup (run from the repository root, e.g. on the login node of the cluster):
#   Rscript setup.R                      # install packages, compile the kernels, run checks
#   Rscript setup.R /path/to/csmGmmCorr  # install csmGmmCorr from a local copy instead of GitHub
# Packages are installed into the first library of .libPaths() (set R_LIBS_USER to choose it).
args <- commandArgs(trailingOnly = TRUE)
repos <- "https://cloud.r-project.org"
lib <- .libPaths()[1]
cran <- c("here", "data.table", "ggplot2", "cowplot", "remotes",
          "qch", "repfdr",                                                   # comparison methods
          "dplyr", "lpSolve", "R.methodsS3", "R.oo", "R.utils", "stringr", "tibble", "tidyr")   # mediation.test
for (p in cran) if (!requireNamespace(p, quietly = TRUE)) install.packages(p, lib = lib, repos = repos)
if (!requireNamespace("adaFilter", quietly = TRUE))
  message("adaFilter is not installed: the comparison uses the built-in implementation in SupportingCode/comparators/adafilter.R")
if (!requireNamespace("mediation.test", quietly = TRUE))
  install.packages("SupportingCode/packages/mediation.test", lib = lib, repos = NULL, type = "source")
if (length(args)) {
  install.packages(args[1], lib = lib, repos = NULL, type = "source")
} else if (!requireNamespace("csmGmmCorr", quietly = TRUE) || utils::packageVersion("csmGmmCorr") < "1.3.0") {
  remotes::install_github(Sys.getenv("CSMGMMCORR_REPO", "zliu20/csmGmmCorr"), lib = lib, upgrade = "never")
}
# compiled kernels
src <- file.path("SupportingCode", "src")
owd <- setwd(src); r <- system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD SHLIB -o",
                                     paste0("reprokernels", .Platform$dynlib.ext), "rd_kernels.cpp glgc_hash.cpp")); setwd(owd)
if (r != 0) stop("compilation of SupportingCode/src failed")

# checks
here::i_am("setup.R")
source(here::here("SupportingCode", "load_all.R"))
load_kernels()
suppressMessages(library(csmGmmCorr))
cat("csmGmmCorr", as.character(packageVersion("csmGmmCorr")), "\n")
for (p in c("qch", "repfdr", "mediation.test", "adaFilter", "data.table", "ggplot2", "cowplot"))
  cat(sprintf("%-15s %s\n", p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "not installed"))
set.seed(1); Z <- matrix(rnorm(3000), 600, 5); s <- runif(600)
for (tg in c("pleiotropy", "replication")) {
  a <- csm_adjust(Z, s, tg, 1L, all = TRUE)$adjusted; b <- rd_majorant(Z, s, tg, 1L)
  rej <- s < 0.3
  ok <- max(abs(a - b)) < 1e-12 && sum(rd_incon_dec(Z, rej, tg, 1L)) == incon_decision(Z, rej, tg, 1L)
  cat(sprintf("kernel check (%s): %s\n", tg, if (ok) "ok" else "FAILED"))
  if (!ok) stop("kernel check failed")
}
h <- .Call("glgc_blake2b8_hex", c("chr1:749963:TAA:T", "chr22:17037192:C:T"))
ok <- identical(h, c("80cbe36f121d036a", "c8e751c24846a6ce"))   # Python hashlib.blake2b(key, digest_size = 8)
cat("BLAKE2b check:", if (ok) "ok" else "FAILED", "\n"); if (!ok) stop("hash check failed")
cat("SETUP DONE\n")
