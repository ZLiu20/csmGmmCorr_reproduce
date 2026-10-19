# Comparison of methods on the lung SCC mediation data, in the spirit of Table 1 of Sun, McCaw
# and Lin (2025): rejections at q = 0.1 (t = 1, 2) for the correlated csmGmm, the csmGmm with
# Sigma = I, DEI, MediationTest and JS-PC, and the top 20 SNP-gene pairs of the correlated
# csmGmm with the number of methods rejecting each (numRej). Run after the fits.
#   Rscript Fig4/mediation_scc_compare.R [all|top3|cis15|cis3]     (default cis15)
here::i_am("Fig4/mediation_scc_compare.R")
source(here::here("SupportingCode", "load_all.R"))
suppressMessages({ library(data.table); library(csmGmmCorr) })
MEDSCC_COLS <- c("alpha_lung", "beta_lung", "alpha_blood", "beta_blood", "alpha_liver", "beta_liver")
snpset <- commandArgs(trailingOnly = TRUE)[1]; if (is.na(snpset)) snpset <- "cis15"   # default: the analysis of the paper
stopifnot(snpset %in% c("all", "top3", "cis15", "cis3"))
cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "4")); q <- 0.1
dd <- here::here("Data", "mediation"); od <- here::here("Fig4", "output", "mediation_scc")
CAP <- as.numeric(Sys.getenv("CAP", "12")); AMU <- as.numeric(Sys.getenv("AMU", "0.5"))
sfx <- paste0(if (AMU == 1.5) "_mu15" else "", if (CAP != 12) paste0("_cap", CAP) else "")   # fits of jobs 10-11 of mediation_scc_fit.R
out <- file.path(od, paste0("compare_", snpset, sfx)); dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- freadgz(file.path(dd, if (startsWith(snpset, "cis")) "mediation_scc_cis.tsv.gz" else "mediation_scc.tsv.gz"))
if (snpset %in% c("top3", "cis3")) d <- d[snp %in% head(fread(file.path(dd, "mediation_scc_snp_rank.tsv"))[order(rank_ILCCO)]$snp, 3)]
d <- d[stats::complete.cases(d[, MEDSCC_COLS, with = FALSE])]
Z <- capz(as.matrix(d[, MEDSCC_COLS, with = FALSE]), CAP); A <- abs(Z)
strong <- (A[, 1] >= 4 & A[, 2] >= 4) | (A[, 3] >= 4 & A[, 4] >= 4) | (A[, 5] >= 4 & A[, 6] >= 4)
weak_alpha <- pmax(A[, 1], A[, 3], A[, 5]) < 2
FITS <- c("Correlated csmGmm" = paste0(snpset, "_main", sfx), "csmGmm, Sigma = I" = paste0(snpset, "_indep", sfx))
COMP <- list(DEI = run_dei_mediation_R1, MediationTest = run_mediation_test_R1, "JS-PC" = run_mediation_safeguard_R1)
for (t in 1:2) {
  rej <- list(); score <- list()
  for (m in names(FITS)) {
    fit <- readRDS(file.path(od, FITS[[m]], "fit.rds"))
    lf <- csm_lfdr(fit, Z, "mediation", t, cores = cores)
    rej[[m]] <- csm_stepup(lf, q)$rejected; score[[m]] <- lf
  }
  for (m in names(COMP)) {
    o <- tryCatch(COMP[[m]](Z, t, q), error = function(e) { message(m, ": ", conditionMessage(e)); NULL })
    if (is.null(o)) next
    rej[[m]] <- as.logical(o$rejected); score[[m]] <- if (!is.null(o$score)) o$score else as.numeric(!o$rejected)
  }
  tab <- rbindlist(lapply(names(rej), function(m) data.table(t = t, method = m, rejected_pairs = sum(rej[[m]]),
    genes = uniqueN(d$gene[rej[[m]]]), rejected_strong_pair = sum(rej[[m]] & strong),
    rejected_all_alpha_lt2 = sum(rej[[m]] & weak_alpha))))
  fwrite(tab, file.path(out, "methods.tsv"), sep = "\t", append = t > 1)
  options(width = 250); print(tab)
  m0 <- names(FITS)[1]
  idx <- which(rej[[m0]]); idx <- head(idx[order(score[[m0]][idx])], 20)
  top <- data.table(rank = seq_along(idx), d[idx, .(snp, locus_gene, gene)], round(Z[idx, , drop = FALSE], 2),
                    lfdr = signif(score[[m0]][idx], 2), strong_pair = strong[idx])
  for (m in names(rej)) top[, (m) := rej[[m]][idx]]
  top[, numRej := rowSums(as.matrix(top[, names(rej), with = FALSE]))]
  fwrite(top, file.path(out, sprintf("top_pairs_t%d.tsv", t)), sep = "\t")
  cat(sprintf("\nt = %d: top 20 SNP-gene pairs of the correlated csmGmm\n", t)); print(top)
}
say("DONE")
