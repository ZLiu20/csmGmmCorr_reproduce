# Second check of the csmGmm comparison against Sun's csmGmm package (upstream HEAD),
#   Rscript Validation/validate2_csmgmm.R          -> Validation/validate2_log.txt
# A. One EM step: from the same parameters, one iteration of csmGmm::symm_fit_ind_EM and one
#    ECM map of csmGmmCorr (Sigma fixed at I, no penalty) must give the same update when no
#    constraint is active (data simulated from the csmGmm, all supports well populated).
#    Same for csmGmm::symm_fit_cor_EM (K = 2, fixed correlation) and fix_sigma = R.
# B. Full fits on the same well-identified data from the same start: log-likelihoods (one
#    common function), lfdr and rejections.
# C. The simulation data of validate_csmgmm.R (pleiotropy, K = 4, t = 3): log-likelihood of
#    both solutions, to see whether they differ in fit or only in poorly identified
#    parameters of nearly empty supports (sizes 3 and 4).
here::i_am("Validation/validate2_csmgmm.R")
if (!requireNamespace("csmGmm", quietly = TRUE)) remotes::install_github("ryanrsun/csmGmm", upgrade = "never")
suppressMessages({ library(csmGmm); library(csmGmmCorr) })
source(here::here("SupportingCode", "load_all.R"))
ns <- asNamespace("csmGmmCorr")
log <- here::here("Validation", "validate2_log.txt"); cat("", file = log)
say <- function(...) { cat(..., "\n", file = log, append = TRUE); cat(..., "\n") }
to_sun <- function(p, sp) {
  K <- sp$K; mu <- list(); pi <- list()
  for (b in 0:(2^K - 1)) {
    g <- which(sp$grp_b == b); g <- g[order(sp$grp_m[g])]
    mu[[b + 1]] <- matrix(t(p$mu[g, , drop = FALSE]), nrow = K)
    pi[[b + 1]] <- p$varpi[g]
  }
  list(mu = mu, pi = pi)
}
from_sun <- function(r, sp, Sigma) {
  G <- length(sp$grp_b); mu <- matrix(0, G, sp$K); v <- numeric(G)
  for (g in seq_len(G)) {
    mu[g, ] <- r$muInfo[[sp$grp_b[g] + 1]][, sp$grp_m[g]]
    v[g] <- r$piInfo[[sp$grp_b[g] + 1]][sp$grp_m[g]]
  }
  list(varpi = v, mu = mu, Sigma = Sigma)
}
loglik <- function(Z, sp, p) ns$.estep(Z, sp, p, block_cells = 5e6, cores = 1L)$F * nrow(Z)
one_ecm <- function(Z, sp, p, ctl) {
  es <- ns$.estep(Z, sp, p, block_cells = 5e6, cores = 1L)
  ns$.ecm_map(es, p, sp, crossprod(Z) / nrow(Z), ctl, nrow(Z))
}
sun_one <- function(Z, s, R = NULL) {           # eps = 9.99 < initial difference 10: exactly one iteration
  if (is.null(R)) symm_fit_ind_EM(Z, s$mu, s$pi, eps = 9.99, checkpoint = FALSE)
  else symm_fit_cor_EM(Z, R, s$mu, s$pi, eps = 9.99, checkpoint = FALSE)
}
truth <- function(sp, null = 0.6) {
  G <- length(sp$grp_b); val <- numeric(G)
  val[sp$grp_size == 1 & sp$grp_m == 1] <- 2; val[sp$grp_size == 1 & sp$grp_m == 2] <- 4
  val[sp$grp_size >= 2] <- 3; val[sp$grp_b == 2^sp$K - 1] <- 5
  v <- ifelse(sp$grp_b == 0, null, (1 - null) / (G - 1))
  list(varpi = v, mu = sp$grp_A * val)
}

# ---- A. one step, K = 4, Sigma = I ----
set.seed(1)
K <- 4L; sp <- csm_spec(K, M = csm_structure(K, 2L, 1L), b_mu = 15)
tr <- truth(sp)
pt <- csm_params(sp, varpi = tr$varpi, mu = tr$mu)
Z <- csm_simulate(20000, pt, sp)$Z
p0 <- list(varpi = 0.7 * tr$varpi + 0.3 / length(tr$varpi), mu = tr$mu * 1.1, Sigma = diag(K))
ctl0 <- csm_control(penalty = 0, weight_penalty = 0, fix_sigma = TRUE)
o1 <- one_ecm(Z, sp, p0, ctl0)
s1 <- sun_one(Z, to_sun(p0, sp))
u <- from_sun(s1, sp, diag(K))
say(sprintf("A1. K = 4, Sigma = I, one step: max |weight difference| %.2e, max |magnitude difference| %.2e (Sun iterations %d)",
            max(abs(o1$varpi - u$varpi)), max(abs(o1$mu - u$mu)), s1$iter))
# ---- A. one step, K = 2, fixed correlation ----
K2 <- 2L; sp2 <- csm_spec(K2, M = csm_structure(K2, 2L, 1L), b_mu = 15)
tr2 <- truth(sp2); R <- matrix(c(1, 0.4, 0.4, 1), 2)
Z2 <- csm_simulate(20000, csm_params(sp2, varpi = tr2$varpi, mu = tr2$mu, Sigma = R), sp2)$Z
q0 <- list(varpi = 0.7 * tr2$varpi + 0.3 / length(tr2$varpi), mu = tr2$mu * 1.1, Sigma = R)
o2 <- one_ecm(Z2, sp2, q0, csm_control(penalty = 0, weight_penalty = 0, fix_sigma = R))
s2 <- sun_one(Z2, to_sun(q0, sp2), R)
u2 <- from_sun(s2, sp2, R)
say(sprintf("A2. K = 2, correlation 0.4 fixed, one step: max |weight difference| %.2e, max |magnitude difference| %.2e",
            max(abs(o2$varpi - u2$varpi)), max(abs(o2$mu - u2$mu))))

# ---- B. full fits from the same start, K = 4, Sigma = I ----
st <- ns$.start_params(Z, sp, "null")
ctl <- csm_control(penalty = 0, weight_penalty = 0, fix_sigma = TRUE, tol_loglik = 1e-9 / nrow(Z), max_iter = 20000)
t1 <- system.time(f <- csm_fit(Z, sp, starts = character(0), init = list(csm_params(sp, varpi = st$varpi, mu = st$mu)), control = ctl))[["elapsed"]]
t2 <- system.time(r <- symm_fit_ind_EM(Z, to_sun(st, sp)$mu, to_sun(st, sp)$pi, eps = 1e-10, checkpoint = FALSE))[["elapsed"]]
ps <- from_sun(r, sp, diag(K))
lf_o <- csm_lfdr(f, Z, "pleiotropy", K - 1L)
rej <- function(l) csm_stepup(pmin(pmax(l, 0), 1), 0.1)$rejected
say(sprintf("B. K = 4, same start: csmGmmCorr %.1f s, log-likelihood %.4f; csmGmm %d iterations, %.1f s, log-likelihood %.4f",
            t1, loglik(Z, sp, f$params), r$iter, t2, loglik(Z, sp, ps)))
say(sprintf("   max |weight difference| %.2e, max |magnitude difference| %.2e, max |lfdr difference| %.2e, rejections %d / %d (both %d)",
            max(abs(f$params$varpi - ps$varpi)), max(abs(f$params$mu - ps$mu)), max(abs(lf_o - r$lfdrResults)),
            sum(rej(lf_o)), sum(rej(r$lfdrResults)), sum(rej(lf_o) & rej(r$lfdrResults))))

# ---- C. the simulation data of validate_csmgmm.R ----
x <- make_dataset("pl", 5L, 4L, 0L, 1L, N = 20000L)$train
Zc <- x$Z
fc <- csm_fit(Zc, sp, starts = c("null", "shared"), control = ctl)
rc <- symm_fit_ind_EM(Zc, to_sun(fc$params, sp)$mu, to_sun(fc$params, sp)$pi, eps = 1e-10, checkpoint = FALSE)
pc <- from_sun(rc, sp, diag(4))
sz <- tabulate(rowSums(x$H != 0) + 1, 5)
say(sprintf("C. simulation data (support sizes %s): csmGmmCorr log-likelihood %.4f; csmGmm from it: %d iterations, log-likelihood %.4f",
            paste(sz, collapse = "/"), loglik(Zc, sp, fc$params), rc$iter, loglik(Zc, sp, pc)))
d <- abs(pc$mu - fc$params$mu); g <- which(d == max(d), arr.ind = TRUE)[1, ]
say(sprintf("   largest magnitude change: support size %d, component %d, weight %.2e (%.3f -> %.3f); smallest csmGmm magnitude %.3f",
            sp$grp_size[g[1]], sp$grp_m[g[1]], fc$params$varpi[g[1]], fc$params$mu[g[1], g[2]], pc$mu[g[1], g[2]],
            min(pc$mu[sp$grp_A])))
say("VALIDATION 2 DONE")
