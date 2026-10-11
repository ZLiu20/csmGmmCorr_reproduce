# Comparison of methods on the LDL mediation application (data set "noshare", K = 6), in the
# spirit of Table 1 of Sun, McCaw and Lin (2025): for the top loci of the correlated csmGmm, the
# z-statistics of the three pathways and the methods that reject the lead variant.
# Methods at q = 0.1, targets t = 1, 2:
#   correlated csmGmm (estimated Sigma), a_mu = 3 and 0.5   (fits noshare_main_amu30, noshare_main)
#   csmGmm with Sigma = I (Sun et al.), a_mu = 3 and 0.5   (fits noshare_indep_amu30, noshare_indep)
#   DEI, MediationTest (minimax) and JS-PC                  (SupportingCode/comparators)
#   Rscript Fig4/mediation_ldl_compare.R            (or bsub < Fig4/run_mediation_ldl_compare.lsf)
# Output in Fig4/output/mediation_ldl/compare/: methods.tsv (rejections, loci, rejections with a
# strong pair, rejections with every |alpha| < 2) and top_loci_t<t>.tsv (top 30 loci).
here::i_am("Fig4/mediation_ldl_compare.R")
source(here::here("SupportingCode", "load_all.R"))
suppressMessages({ library(data.table); library(csmGmmCorr) })
cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "4")); q <- 0.1
od <- here::here("Fig4", "output", "mediation_ldl"); out <- file.path(od, "compare")
dir.create(out, showWarnings = FALSE)
d <- freadgz(here::here("Data", "mediation", "mediation_ldl_liver_noshare.tsv.gz"))
d <- d[stats::complete.cases(d[, MEDLDL_COLS, with = FALSE])]
Z <- capz(as.matrix(d[, MEDLDL_COLS, with = FALSE]), 12)
A <- abs(Z)
strong <- (A[, 1] >= 4 & A[, 2] >= 4) | (A[, 3] >= 4 & A[, 4] >= 4) | (A[, 5] >= 4 & A[, 6] >= 4)
weak_alpha <- pmax(A[, 1], A[, 3], A[, 5]) < 2
chr <- as.character(d$chr); pos <- d$pos_b38
FITS <- c("Correlated csmGmm (a_mu 3)" = "noshare_main_amu30", "csmGmm, Sigma = I (a_mu 3)" = "noshare_indep_amu30",
          "Correlated csmGmm (a_mu 0.5)" = "noshare_main", "csmGmm, Sigma = I (a_mu 0.5)" = "noshare_indep")
COMP <- list(DEI = run_dei_mediation_R1, MediationTest = run_mediation_test_R1, "JS-PC" = run_mediation_safeguard_R1)
loci_n <- function(sel, score) if (any(sel)) length(clump_leads(chr, pos, score, sel)) else 0L
for (t in 1:2) {
  rej <- list(); score <- list()
  for (m in names(FITS)) {
    fit <- readRDS(file.path(od, FITS[[m]], "fit.rds"))
    lf <- csm_lfdr(fit, Z, "mediation", t, cores = cores)
    rej[[m]] <- csm_stepup(lf, q)$rejected; score[[m]] <- lf
    say("t = %d, %s: %d", t, m, sum(rej[[m]]))
  }
  for (m in names(COMP)) {
    o <- tryCatch(COMP[[m]](Z, t, q), error = function(e) { message(m, ": ", conditionMessage(e)); NULL })
    if (is.null(o)) next
    rej[[m]] <- as.logical(o$rejected); score[[m]] <- if (!is.null(o$score)) o$score else as.numeric(!o$rejected)
    say("t = %d, %s: %d", t, m, sum(rej[[m]]))
  }
  tab <- rbindlist(lapply(names(rej), function(m) data.table(t = t, method = m, rejected = sum(rej[[m]]),
    loci = loci_n(rej[[m]], score[[m]]), rejected_strong_pair = sum(rej[[m]] & strong),
    rejected_all_alpha_lt2 = sum(rej[[m]] & weak_alpha))))
  fwrite(tab, file.path(out, "methods.tsv"), sep = "\t", append = t > 1)
  options(width = 250); print(tab)
  # top loci of the correlated csmGmm (a_mu = 3): lead variants, with the methods that reject them
  m0 <- names(FITS)[1]
  lead <- clump_leads(chr, pos, score[[m0]], rej[[m0]])
  lead <- head(lead[order(score[[m0]][lead])], 30)
  top <- data.table(rank = seq_along(lead), chr = d$chr[lead], pos_b38 = pos[lead],
                    gene1 = d$gene1_name[lead], gene2 = d$gene2_name[lead], gene3 = d$gene3_name[lead],
                    round(Z[lead, , drop = FALSE], 1), lfdr = signif(score[[m0]][lead], 2), strong_pair = strong[lead])
  for (m in names(rej)) top[, (m) := rej[[m]][lead]]
  top[, numRej := rowSums(as.matrix(top[, names(rej), with = FALSE]))]
  fwrite(top, file.path(out, sprintf("top_loci_t%d.tsv", t)), sep = "\t")
  SHORT <- c("Correlated csmGmm (a_mu 3)" = "CcsmGmm3", "csmGmm, Sigma = I (a_mu 3)" = "csmGmmI3",
             "Correlated csmGmm (a_mu 0.5)" = "CcsmGmm05", "csmGmm, Sigma = I (a_mu 0.5)" = "csmGmmI05",
             DEI = "DEI", MediationTest = "MedTest", "JS-PC" = "JSPC")
  short <- unname(SHORT[names(rej)])        # labels follow the methods that ran (a failed comparator is absent)
  pr <- copy(top); setnames(pr, names(rej), short)
  for (s in short) pr[, (s) := fifelse(get(s), "x", "")]
  cat(sprintf("\nt = %d: top 30 loci of the correlated csmGmm (a_mu = 3); x = lead variant rejected\n", t))
  print(pr[, c("rank", "chr", "gene1", "gene2", "gene3", MEDLDL_COLS, "lfdr", "strong_pair", short, "numRej"), with = FALSE], nrows = 30)
}
say("DONE")
