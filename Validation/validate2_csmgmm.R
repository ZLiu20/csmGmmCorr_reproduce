here::i_am("Validation/validate2_csmgmm.R")
suppressPackageStartupMessages(library(csmGmmCorr))
source(here::here("SupportingCode", "load_all.R"))

log_file <- here::here("Validation", "validate2_log.txt")
cat("", file = log_file)
say <- function(...) {
  line <- paste0(..., collapse = "")
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

N <- as.integer(Sys.getenv("VAL_N", "20000"))
if (!is.finite(N) || N < 500L) stop("VAL_N must be at least 500")
if (!requireNamespace("csmGmm", quietly = TRUE)) {
  say("EXTERNAL_CSMGMM_SKIPPED package_not_installed")
  quit(save = "no", status = 0L)
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
  list(varpi = weights, mu = spec$grp_A * level)
}

one_ours <- function(z, spec, params, fixed_sigma) {
  control <- csm_control(
    max_iter = 2L, accelerate = FALSE, penalty = 0,
    weight_penalty = 0, fix_sigma = fixed_sigma
  )
  es <- csmGmmCorr:::.estep(z, spec, params, cores = 1L)
  csmGmmCorr:::.ecm_map(
    es, params, spec, crossprod(z) / nrow(z), control, nrow(z)
  )
}

one_sun <- function(z, params, spec, R = NULL) {
  initial <- as_sun_lists(params, spec)
  if (is.null(R)) {
    csmGmm::symm_fit_ind_EM(
      z, initial$muList, initial$piList,
      eps = 9.99, checkpoint = FALSE
    )
  } else {
    csmGmm::symm_fit_cor_EM(
      z, R, initial$muList, initial$piList,
      eps = 9.99, checkpoint = FALSE
    )
  }
}

from_sun <- function(result, spec, Sigma) {
  csm_params(spec, pi = result$piInfo, mu = result$muInfo, Sigma = Sigma)
}

compare_one_step <- function(K, R, seed) {
  set.seed(seed)
  spec <- csm_spec(K, M = csm_structure(K, 2L, 1L), S = K, b_mu = 15)
  truth <- make_truth(spec)
  true_par <- csm_params(spec, varpi = truth$varpi, mu = truth$mu, Sigma = R)
  z <- csm_simulate(N, true_par, spec)$Z
  G <- length(truth$varpi)
  initial <- csm_params(
    spec,
    varpi = 0.7 * truth$varpi + 0.3 / G,
    mu = 1.1 * truth$mu,
    Sigma = R
  )
  ours <- one_ours(z, spec, initial, R)
  upstream <- one_sun(z, initial, spec, if (K == 4L) NULL else R)
  reference <- from_sun(upstream, spec, R)
  c(weight = max(abs(ours$varpi - reference$varpi)),
    magnitude = max(abs(ours$mu - reference$mu)),
    iterations = upstream$iter)
}

independent <- compare_one_step(4L, diag(4L), 7401)
correlated <- compare_one_step(2L, matrix(c(1, 0.4, 0.4, 1), 2L), 7402)
say(sprintf("one step K=4: weight %.3e; magnitude %.3e; upstream iterations %d",
            independent["weight"], independent["magnitude"],
            independent["iterations"]))
say(sprintf("one step K=2, rho=.4: weight %.3e; magnitude %.3e; upstream iterations %d",
            correlated["weight"], correlated["magnitude"],
            correlated["iterations"]))
stopifnot(independent["weight"] < 1e-10,
          independent["magnitude"] < 1e-10,
          correlated["weight"] < 1e-10,
          correlated["magnitude"] < 1e-10)

simulation_n <- max(500L, min(N, 2000L))
generated <- make_dataset("pl", 5L, 4L, 0L, 1L, N = simulation_n)$train
z <- generated$Z
spec <- csm_spec(4L, M = csm_structure(4L, 2L, 1L), S = 4L,
                 a_mu = 1e-8, b_mu = 15)
start <- csmGmmCorr:::.start_params(z, spec, "null")
maps <- as.integer(Sys.getenv("VAL_MAPS", "8"))
if (!is.finite(maps) || maps < 1L) stop("VAL_MAPS must be positive")
ours_par <- start
upstream_par <- start
for (iteration in seq_len(maps)) {
  ours_par <- one_ours(z, spec, ours_par, diag(4L))
  upstream <- one_sun(z, upstream_par, spec)
  upstream_par <- from_sun(upstream, spec, diag(4L))
}
upstream_lfdr <- csm_lfdr(upstream_par, z, "pleiotropy", 3L, spec = spec)
ours_lfdr <- csm_lfdr(ours_par, z, "pleiotropy", 3L, spec = spec)
agreement <- suppressWarnings(cor(ours_lfdr, upstream_lfdr))
say(sprintf(paste0("generated data: maps=%d; corr(lfdr)=%.6f; ",
                   "max weight delta=%.3e; max magnitude delta=%.3e"),
            maps, agreement,
            max(abs(ours_par$varpi - upstream_par$varpi)),
            max(abs(ours_par$mu - upstream_par$mu))))
stopifnot(is.finite(agreement), agreement > 0.99,
          max(abs(ours_par$varpi - upstream_par$varpi)) < 0.05,
          max(abs(ours_par$mu - upstream_par$mu)) < 1)

say("VALIDATION_2_DONE")
