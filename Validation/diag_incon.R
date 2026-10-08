# Diagnosis of the incongruous results of the correlated csmGmm :
#   Rscript Validation/diag_incon.R            -> Validation/diag_incon_log.txt
# For one replicate of three designs with many incongruous results (K = 4, t = 1):
#  1. our count (R/incongruity.R) and csmGmm::check_incongruous (Sun's code) on the same scores;
#  2. the flagged features: size of the score differences, truncated statistics, sign patterns;
#  3. the same data refitted with Sigma = I (independence), to separate the effect of the
#     correlation from that of the model structure and the target (t < K - 1);
#  4. the target t = K - 1, for which Sun et al.'s ordering result applies.
here::i_am("Validation/diag_incon.R")
if (!requireNamespace("csmGmm", quietly = TRUE)) remotes::install_github("ryanrsun/csmGmm", upgrade = "never")
suppressMessages({ library(csmGmm); library(csmGmmCorr) })
source(here::here("SupportingCode", "load_all.R"))
log <- here::here("Validation", "diag_incon_log.txt"); cat("", file = log)
say <- function(...) { cat(..., "\n", file = log, append = TRUE); cat(..., "\n") }
N <- as.integer(Sys.getenv("SIM_N", "100000")); cores <- as.integer(Sys.getenv("CORES", "4"))

# flagged features (Sun's definition) with the largest score difference to a dominated feature
incon_detail <- function(Z, s, target, t, max_n = 5000L) {
  keep <- which(is.finite(s) & s < 0.99)
  keep <- keep[order(s[keep])][seq_len(min(length(keep), max_n))]
  rp <- .incon_repr(Z[keep, , drop = FALSE], target, t); W <- rp$w; K <- ncol(W); sc <- s[keep]
  out <- NULL
  for (gr in unique(rp$group)) {
    m <- which(rp$group == gr); if (length(m) < 2L) next
    m <- m[do.call(order, c(list(sc[m]), lapply(seq_len(K), function(k) -W[m, k])))]
    for (i in 2:length(m)) {
      prev <- m[seq_len(i - 1L)]; wi <- rep(W[m[i], ], each = length(prev))
      dom <- rowSums(W[prev, , drop = FALSE] <= wi) == K & rowSums(W[prev, , drop = FALSE] < wi) > 0
      if (any(dom)) {
        j <- prev[dom][which.min(sc[prev[dom]])]
        out <- rbind(out, data.frame(i = keep[m[i]], j = keep[j], s_i = sc[m[i]], s_j = sc[j]))
      }
    }
  }
  out
}
fitB <- function(Z, K, fix = FALSE) {
  sp <- csm_spec(K, M = csm_structure(K, 2L, 1L), b_mu = 15)
  csm_fit(Z, sp, starts = c("null", "shared"), control = csm_control(cores = cores, max_iter = 1500, fix_sigma = fix))
}
cases <- list(list(scen = "pl", d = 7L, target = "pleiotropy", name = "pleiotropy tau1 = 0.02"),
              list(scen = "rp", d = 9L, target = "replication", name = "replication tau1 = 0.04"),
              list(scen = "med", d = 6L, target = "mediation", name = "mediation tau1 = 0.02"))
K <- 4L; t <- 1L
for (cs in cases) {
  dat <- make_dataset(cs$scen, cs$d, K, if (cs$scen == "med") t else 0L, 1L, N)$train
  Z <- dat$Z
  say(sprintf("==== %s, K = %d, t = %d, replicate 1 ====", cs$name, K, t))
  f <- fitB(Z, K)
  say(sprintf("fit: converged %s; Sigma off-diagonal %s", f$converged,
              paste(sprintf("%.3f", f$params$Sigma[upper.tri(f$params$Sigma)]), collapse = " ")))
  lf <- csm_lfdr(f, Z, cs$target, t, cores = cores)
  rej <- csm_stepup(lf, 0.1)$rejected
  ours <- incon_score(Z, lf, cs$target, t)
  keep <- which(lf < 0.99); keep <- keep[order(lf[keep])][seq_len(min(length(keep), 5000L))]
  sun <- tryCatch(length(check_incongruous(Z[keep, ], lf[keep])), error = function(e) NA)
  say(sprintf("1. incongruous: ours %d; csmGmm::check_incongruous on the same 5000 scores %s (orthant order; differs for replication by design)",
              ours, as.character(sun)))
  dd <- incon_detail(Z, lf, cs$target, t)
  if (!is.null(dd) && nrow(dd)) {
    gap <- dd$s_i - dd$s_j
    trunc <- rowSums(abs(Z[dd$i, , drop = FALSE]) >= 8.1 - 1e-9) > 0
    say(sprintf("2. flagged %d: score difference quantiles (50/90/100%%) %.2e / %.2e / %.2e; below 1e-6: %d; lfdr of flagged median %.3f; rejected %d; truncated |z| = 8.1: %d",
                nrow(dd), quantile(gap, .5), quantile(gap, .9), max(gap), sum(gap < 1e-6), median(dd$s_i), sum(rej[dd$i]), sum(trunc)))
    sz <- rowSums(dat$H[dd$i, , drop = FALSE] != 0)
    say(sprintf("   true support sizes of flagged features: %s", paste(names(table(sz)), table(sz), sep = ":", collapse = " ")))
    ex <- head(dd[order(-(dd$s_i - dd$s_j)), ], 5)
    for (r in seq_len(nrow(ex))) say(sprintf("   e.g. z_i = (%s), lfdr %.4f  dominates  z_j = (%s), lfdr %.4f",
      paste(sprintf("%.2f", Z[ex$i[r], ]), collapse = ", "), ex$s_i[r], paste(sprintf("%.2f", Z[ex$j[r], ]), collapse = ", "), ex$s_j[r]))
  }
  fI <- fitB(Z, K, fix = TRUE)
  lfI <- csm_lfdr(fI, Z, cs$target, t, cores = cores)
  say(sprintf("3. Sigma fixed at I: incongruous %d (rejections %d vs %d with Sigma estimated)",
              incon_score(Z, lfI, cs$target, t), sum(csm_stepup(lfI, 0.1)$rejected), sum(rej)))
  tg2 <- if (cs$target == "mediation") "mediation" else cs$target
  lf3 <- csm_lfdr(f, Z, tg2, K - 1L, cores = cores); lf3I <- csm_lfdr(fI, Z, tg2, K - 1L, cores = cores)
  say(sprintf("4. target t = K - 1 = %d: incongruous %d (Sigma estimated), %d (Sigma = I)",
              K - 1L, incon_score(Z, lf3, tg2, K - 1L), incon_score(Z, lf3I, tg2, K - 1L)))
}
say("DIAGNOSIS DONE")
