here::i_am("Validation/validate_csmgmm.R")
suppressPackageStartupMessages(library(csmGmmCorr))
source(here::here("SupportingCode", "load_all.R"))

log_file <- here::here("Validation", "validate_csmgmm_log.txt")
cat("", file = log_file)
say <- function(...) {
  line <- paste0(..., collapse = "")
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

N <- as.integer(Sys.getenv("VAL_N", "20000"))
if (!is.finite(N) || N < 500L) stop("VAL_N must be at least 500")
fit_tolerance <- as.numeric(Sys.getenv("VAL_TOL", "1e-4"))
if (!is.finite(fit_tolerance) || fit_tolerance <= 0) {
  stop("VAL_TOL must be positive")
}

make_truth <- function(spec, null_mass = 0.6) {
  G <- length(spec$grp_b)
  weights <- rep((1 - null_mass) / (G - 1L), G)
  weights[spec$grp_b == 0L] <- null_mass
  level <- rep(3, G)
  level[spec$grp_b == 0L] <- 0
  level[spec$grp_size == 1L & spec$grp_m == 1L] <- 2
  level[spec$grp_size == 1L & spec$grp_m == 2L] <- 4
  level[spec$grp_b == 2^spec$K - 1L] <- 5
  csm_params(spec, varpi = weights, mu = spec$grp_A * level,
             Sigma = diag(spec$K))
}

from_sun <- function(result, spec, Sigma) {
  csm_params(spec, pi = result$piInfo, mu = result$muInfo, Sigma = Sigma)
}

fit_from <- function(z, spec, initial, fixed_sigma) {
  csm_fit(
    z, spec, starts = character(), init = list(initial),
    control = csm_control(
      max_iter = 5000L, tol_loglik = fit_tolerance, accelerate = TRUE,
      penalty = 0, weight_penalty = 0, fix_sigma = fixed_sigma
    )
  )
}

if (!requireNamespace("csmGmm", quietly = TRUE)) {
  say("EXTERNAL_CSMGMM_SKIPPED package_not_installed")
  quit(save = "no", status = 0L)
}
say("csmGmm ", as.character(utils::packageVersion("csmGmm")),
    "; csmGmmCorr ", as.character(utils::packageVersion("csmGmmCorr")))

set.seed(7301)
K <- 4L
spec <- csm_spec(K, M = csm_structure(K, 2L, 1L), S = K, b_mu = 15)
truth <- make_truth(spec)
z <- csm_simulate(N, truth, spec)$Z
ours <- csm_fit(
  z, spec, starts = c("null", "shared"),
  control = csm_control(
    max_iter = 5000L, tol_loglik = fit_tolerance, accelerate = TRUE,
    penalty = 0, weight_penalty = 0, fix_sigma = TRUE
  )
)
stopifnot(isTRUE(ours$converged))
sun_start <- as_sun_lists(ours)
fixed <- csmGmm::symm_fit_ind_EM(
  z, sun_start$muList, sun_start$piList,
  eps = fit_tolerance, checkpoint = FALSE
)
fixed_par <- from_sun(fixed, spec, diag(K))
fixed_delta <- max(abs(c(ours$params$varpi - fixed_par$varpi,
                         ours$params$mu - fixed_par$mu)))
fixed_lfdr <- csm_lfdr(ours, z, "pleiotropy", K - 1L)
lfdr_delta <- max(abs(fixed_lfdr - fixed$lfdrResults))
say(sprintf("fixed point: upstream iterations=%d; parameter delta=%.3e; lfdr delta=%.3e",
            fixed$iter, fixed_delta, lfdr_delta))
stopifnot(fixed_delta < 0.03, lfdr_delta < 0.03)

initial <- csmGmmCorr:::.start_params(z, spec, "null")
same_ours <- fit_from(z, spec, initial, TRUE)
initial_sun <- as_sun_lists(initial, spec)
same_sun <- csmGmm::symm_fit_ind_EM(
  z, initial_sun$muList, initial_sun$piList,
  eps = fit_tolerance, checkpoint = FALSE
)
same_par <- from_sun(same_sun, spec, diag(K))
same_lfdr <- csm_lfdr(same_ours, z, "pleiotropy", K - 1L)
lfdr_correlation <- suppressWarnings(cor(same_lfdr, same_sun$lfdrResults))
weight_delta <- max(abs(same_ours$params$varpi - same_par$varpi))
magnitude_delta <- max(abs(same_ours$params$mu - same_par$mu))
say(sprintf(paste0("same start: corr(lfdr)=%.6f; max weight delta=%.3e; ",
                   "max magnitude delta=%.3e"),
            lfdr_correlation, weight_delta, magnitude_delta))
stopifnot(is.finite(lfdr_correlation), lfdr_correlation > 0.995,
          weight_delta < 0.03, magnitude_delta < 0.15)

set.seed(7302)
K2 <- 2L
spec2 <- csm_spec(K2, M = csm_structure(K2, 2L, 1L), S = K2, b_mu = 15)
R2 <- matrix(c(1, 0.4, 0.4, 1), 2L)
truth2 <- make_truth(spec2)
truth2 <- csm_params(spec2, varpi = truth2$varpi, mu = truth2$mu, Sigma = R2)
z2 <- csm_simulate(N, truth2, spec2)$Z
ours2 <- fit_from(z2, spec2, truth2, R2)
sun2_start <- as_sun_lists(ours2)
fixed2 <- csmGmm::symm_fit_cor_EM(
  z2, R2, sun2_start$muList, sun2_start$piList,
  eps = fit_tolerance, checkpoint = FALSE
)
fixed2_par <- from_sun(fixed2, spec2, R2)
delta2 <- max(abs(c(ours2$params$varpi - fixed2_par$varpi,
                    ours2$params$mu - fixed2_par$mu)))
lfdr2 <- csm_lfdr(ours2, z2, "pleiotropy", 1L)
lfdr2_delta <- max(abs(lfdr2 - fixed2$lfdrResults))
say(sprintf("fixed correlation: upstream iterations=%d; parameter delta=%.3e; lfdr delta=%.3e",
            fixed2$iter, delta2, lfdr2_delta))
stopifnot(delta2 < 0.03, lfdr2_delta < 0.03)

say("VALIDATION_DONE")
