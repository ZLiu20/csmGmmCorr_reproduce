# Check that the comparison methods are the estimators of Sun, McCaw and Lin's csmGmm
# package (github ryanrsun/csmGmm, version 0.5.0, commit 4a03af5):
#   "csmGmm"     (csmGmmCorr, Sigma = I, ML, structure B)  vs  symm_fit_ind_EM   (K = 4)
#   "csmGmm_cor" (csmGmmCorr, Sigma = plug-in R, ML)       vs  symm_fit_cor_EM   (K = 2, the
#                case it is written for; corMat = cor(testStats) as in its example) 
#   Rscript Validation/validate_csmgmm.R      (installs csmGmm from GitHub, ryanrsun/csmGmm, if needed)
# Writes Validation/validate_csmgmm_log.txt. The target is the one csmGmm supports: all K non-null
# (pleiotropy with t = K - 1).
here::i_am("Validation/validate_csmgmm.R")
if (!requireNamespace("csmGmm", quietly = TRUE)) remotes::install_github("ryanrsun/csmGmm@4a03af5", upgrade = "never")
suppressMessages({ library(csmGmm); library(csmGmmCorr) })
source(here::here("SupportingCode", "load_all.R"))
log <- here::here("Validation", "validate_csmgmm_log.txt"); cat("", file = log)
say <- function(...) { cat(..., "\n", file = log, append = TRUE); cat(..., "\n") }
say(sprintf("csmGmm %s; csmGmmCorr %s", packageVersion("csmGmm"), packageVersion("csmGmmCorr")))

N <- as.integer(Sys.getenv("VAL_N", "20000")); K <- 4L
x <- make_dataset("pl", 5L, K, 0L, 1L, N = N)$train
Z <- x$Z; H <- x$H
sp <- csm_spec(K, M = csm_structure(K, 2L, 1L), b_mu = 15)
ctl <- csm_control(penalty = 0, weight_penalty = 0, fix_sigma = TRUE, max_iter = 5000)
tm1 <- system.time(f <- csm_fit(Z, sp, starts = c("null", "shared"), control = ctl))[["elapsed"]]
say(sprintf("csmGmmCorr (Sigma = I, ML, structure B): N = %d, K = %d, %.1f s, converged %s, Sigma = I: %s",
            N, K, tm1, f$converged, isTRUE(all.equal(f$params$Sigma, diag(K)))))
lf_ours <- csm_lfdr(f, Z, "pleiotropy", K - 1L)

# parameters of our fit as csmGmm lists (support label b: coordinate 1 is the most significant bit)
to_sun <- function(p, sp) {
  mu <- list(); pi <- list()
  K <- sp$K
  for (b in 0:(2^K - 1)) {
    g <- which(sp$grp_b == b); g <- g[order(sp$grp_m[g])]
    mu[[b + 1]] <- matrix(t(p$mu[g, , drop = FALSE]), nrow = K)
    pi[[b + 1]] <- p$varpi[g]
  }
  list(mu = mu, pi = pi)
}
# 1. Fixed point: csmGmm's EM started at our estimate should not move.
s0 <- to_sun(f$params, sp)
tm2 <- system.time(r1 <- symm_fit_ind_EM(Z, s0$mu, s0$pi, eps = 1e-10, checkpoint = FALSE))[["elapsed"]]
d1 <- max(abs(c(unlist(r1$muInfo), unlist(r1$piInfo)) - c(unlist(s0$mu), unlist(s0$pi))))
say(sprintf("1. csmGmm EM started at our estimate: %d iterations, %.1f s, max parameter change %.2e, max |lfdr difference| %.2e",
            r1$iter, tm2, d1, max(abs(r1$lfdrResults - lf_ours))))
# 2. Both fitted from the same "null" starting value.
st <- csmGmmCorr:::.start_params(Z, sp, "null"); s1 <- to_sun(st, sp)
tm3 <- system.time(r2 <- symm_fit_ind_EM(Z, s1$mu, s1$pi, eps = 1e-8, checkpoint = FALSE))[["elapsed"]]
alt <- rowSums(H != 0) == K
rej <- function(l) csm_stepup(pmin(pmax(l, 0), 1), 0.1)$rejected
ro <- rej(lf_ours); rs <- rej(r2$lfdrResults)
say(sprintf("2. csmGmm EM from the same start: %d iterations, %.1f s; max |lfdr difference| %.2e, correlation %.6f",
            r2$iter, tm3, max(abs(r2$lfdrResults - lf_ours)), cor(r2$lfdrResults, lf_ours)))
say(sprintf("   rejections at q = 0.1: ours %d (FDP %.3f), csmGmm %d (FDP %.3f), both %d",
            sum(ro), if (sum(ro)) mean(!alt[ro]) else 0, sum(rs), if (sum(rs)) mean(!alt[rs]) else 0, sum(ro & rs)))
# 3. Plug-in correlation, K = 2 (symm_fit_cor_EM is written for K = 2).
Z2 <- Z[, 1:2]
sp2 <- csm_spec(2L, M = csm_structure(2L, 2L, 1L), b_mu = 15)
R2 <- csmGmmCorr:::.project_cor(cor(Z2), sp2$kappa)
f2 <- csm_fit(Z2, sp2, starts = c("null", "shared"),
              control = csm_control(penalty = 0, weight_penalty = 0, fix_sigma = R2, max_iter = 5000))
K <- 2L; s2 <- to_sun(f2$params, sp2)
tm4 <- system.time(r3 <- symm_fit_cor_EM(Z2, R2, s2$mu, s2$pi, eps = 1e-10, checkpoint = FALSE))[["elapsed"]]
d3 <- max(abs(c(unlist(r3$muInfo), unlist(r3$piInfo)) - c(unlist(s2$mu), unlist(s2$pi))))
lf2 <- csm_lfdr(f2, Z2, "pleiotropy", 1L)
say(sprintf("3. K = 2, plug-in correlation %.3f: csmGmm::symm_fit_cor_EM started at our estimate: %d iterations, %.1f s, max parameter change %.2e, max |lfdr difference| %.2e",
            R2[1, 2], r3$iter, tm4, d3, max(abs(r3$lfdrResults - lf2))))
say("VALIDATION DONE")
