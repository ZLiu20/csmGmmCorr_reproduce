# All targets t = 1, ..., 5 (K = 6, G = 3) for every LDL mediation fit, from the saved fits
# (no refitting). For t >= G = 3 the mediation null equals the pleiotropy null "at most t nonzero"
# (main-text Section 2.1); t = 1, 2 are the mediation-specific targets.
#   Rscript Fig4/mediation_ldl_all_t.R <index>    one fit folder (sorted list), writes all_t.tsv into it
#   bsub < Fig4/run_mediation_ldl_all_t.lsf       all fit folders
#   Rscript Fig4/mediation_ldl_all_t.R summary    combined table
# Columns: rejections at q = 0.1 (raw lfdr and monotone-adjusted); strong_target_* = variants with
# a strong pair (|alpha_g| >= 4 and |beta_g| >= 4 for some g) and more than t coordinates with
# |z| >= 4, i.e. clear evidence for the target, in the data and among the rejections;
# all_alpha_lt2 = rejections with every |alpha_g| < 2; loci = 500-kb loci of the rejections and
# loci_with_strong = those containing a rejected strong-target variant.
here::i_am("Fig4/mediation_ldl_all_t.R")
source(here::here("SupportingCode", "load_all.R"))
suppressMessages({ library(data.table); library(csmGmmCorr) })
od <- here::here("Fig4", "output", "mediation_ldl")
dirs <- sort(medldl_fit_dirs(od)); dirs <- dirs[file.exists(file.path(dirs, "fit.rds"))]
arg <- commandArgs(trailingOnly = TRUE)[1]
if (identical(arg, "summary")) {
  x <- rbindlist(lapply(file.path(dirs, "all_t.tsv"), function(f) if (file.exists(f)) fread(f)), fill = TRUE)
  options(width = 200); print(x, nrows = 500)
  quit(save = "no")
}
loci_of <- function(chr, pos, score, win = 5e5) {
  lead <- clump_leads(chr, pos, score, rep(TRUE, length(pos)), win)
  id <- integer(length(pos))
  for (c0 in unique(chr)) {
    i <- which(chr == c0); l <- lead[chr[lead] == c0]
    id[i] <- l[vapply(pos[i], function(p) which.min(abs(pos[l] - p)), 1L)]
  }
  id
}
if (as.integer(arg) > length(dirs)) quit(save = "no")
v <- dirs[as.integer(arg)]; name <- basename(v)
set <- sub("_(main|indep|robust)(_amu[0-9]+)?$", "", name)
d <- freadgz(here::here("Data", "mediation", sprintf("mediation_ldl_liver_%s.tsv.gz", sub("_dedup$", "", set))))
if (grepl("_dedup$", set)) d <- medldl_dedup(d)
cols <- MEDLDL_COLS
d <- d[stats::complete.cases(d[, cols, with = FALSE])]
Z <- capz(as.matrix(d[, cols, with = FALSE]), 12)
fit <- readRDS(file.path(v, "fit.rds"))
cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "4"))
say("%s (data set %s): %d variants", name, set, nrow(Z))
A <- abs(Z)
pair <- (A[, 1] >= 4 & A[, 2] >= 4) | (A[, 3] >= 4 & A[, 4] >= 4) | (A[, 5] >= 4 & A[, 6] >= 4)
n4 <- rowSums(A >= 4); amax <- pmax(A[, 1], A[, 3], A[, 5])
chr <- as.character(d$chr); pos <- d$pos_b38
rows <- list()
for (t in 1:5) {
  lf <- csm_lfdr(fit, Z, "mediation", t, cores = cores)
  su <- csm_stepup(lf, 0.1); ad <- csm_adjust(Z, lf, "mediation", t, 0.1)
  r <- su$rejected; st <- pair & n4 > t
  loc <- if (any(r)) loci_of(chr[r], pos[r], lf[r]) else integer(0)
  rows[[t]] <- data.table(fit = name, t = t, n = nrow(Z), R_raw = sum(r), R_adj = sum(ad$rejected),
    strong_target_data = sum(st),
    strong_target_loci_data = if (any(st)) length(clump_leads(chr[st], pos[st], -n4[st], rep(TRUE, sum(st)))) else 0L,
    strong_target_rej = sum(r & st), all_alpha_lt2 = sum(r & amax < 2),
    loci = uniqueN(loc), loci_with_strong = uniqueN(loc[st[r]]))
  print(rows[[t]])
}
fwrite(rbindlist(rows), file.path(v, "all_t.tsv"), sep = "\t")
say("DONE %s", name)
