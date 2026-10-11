args <- commandArgs(trailingOnly = TRUE)
package_source <- if (length(args)) args[[1L]] else Sys.getenv("CSMGMMCORR_SRC", "")
dir.create("logs", showWarnings = FALSE)

cran <- c("here", "data.table", "ggplot2", "cowplot", "remotes", "qch", "repfdr",
          "RSQLite", "dplyr", "lpSolve", "R.methodsS3", "R.oo", "R.utils",
          "stringr", "tibble", "tidyr")
missing <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")

if (!requireNamespace("adaFilter", quietly = TRUE)) {
  message("adaFilter is not installed; the bundled compatible implementation will be used")
}

if (!requireNamespace("mediation.test", quietly = TRUE)) {
  install.packages(file.path("SupportingCode", "packages", "mediation.test"),
                   repos = NULL, type = "source")
}

need_method <- !requireNamespace("csmGmmCorr", quietly = TRUE) ||
  packageVersion("csmGmmCorr") < package_version("1.3.0")
if (need_method) {
  if (nzchar(package_source)) {
    install.packages(package_source, repos = NULL, type = "source")
  } else {
    remotes::install_github("ZLiu20/csmGmmCorr")
  }
}
if (packageVersion("csmGmmCorr") != package_version("1.3.0")) {
  warning("The paper used csmGmmCorr 1.3.0; installed version is ", packageVersion("csmGmmCorr"))
}

if (!requireNamespace("csmGmm", quietly = TRUE)) {
  try(remotes::install_github("ryanrsun/csmGmm@4a03af5"), silent = TRUE)
}

writeLines(c(capture.output(sessionInfo()), "", capture.output(installed.packages()[, c("Package", "Version")])),
           file.path("logs", "versions.txt"))

src <- file.path("SupportingCode", "src")
old <- setwd(src)
on.exit(setwd(old), add = TRUE)
objects <- c("rd_kernels.cpp", "glgc_hash.cpp")
status <- system2(file.path(R.home("bin"), "R"), c("CMD", "SHLIB", objects, "-o", paste0("reprokernels", .Platform$dynlib.ext)))
if (status != 0L) stop("Native helper compilation failed")
setwd(old)

source(file.path("SupportingCode", "load_all.R"))
load_kernels(getwd())
set.seed(11)
z <- matrix(rnorm(80), 20, 4)
lfdr <- runif(20)
for (target in c("pleiotropy", "replication")) {
  reference <- csm_adjust(z, lfdr, target, 1L, 0.1)$adjusted
  native <- .Call("rd_majorant", z, lfdr, target, 1L, 64L)
  stopifnot(max(abs(reference - native)) < 1e-10)
}
cat("SETUP_OK\n")
