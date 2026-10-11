here::i_am("Validation/diag_incon.R")
suppressPackageStartupMessages(library(csmGmmCorr))
source(here::here("SupportingCode", "load_all.R"))

log_file <- here::here("Validation", "diag_incon_log.txt")
cat("", file = log_file)
say <- function(...) {
  line <- paste0(..., collapse = "")
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

N <- as.integer(Sys.getenv("SIM_N", "100000"))
cores <- as.integer(Sys.getenv("CORES", "1"))
max_iter <- as.integer(Sys.getenv("DIAG_MAX_ITER", "1500"))
if (!is.finite(N) || N < 500L) stop("SIM_N must be at least 500")
if (!is.finite(cores) || cores < 1L) stop("CORES must be positive")
if (!is.finite(max_iter) || max_iter < 2L) stop("DIAG_MAX_ITER must be at least two")

fit_case <- function(z, K, fixed = FALSE) {
  spec <- csm_spec(K, M = csm_structure(K, 2L, 1L), S = K, b_mu = 15)
  csm_fit(
    z, spec, starts = c("null", "shared"),
    control = csm_control(
      cores = cores, max_iter = max_iter,
      fix_sigma = fixed
    )
  )
}

cases <- list(
  list(scen = "pl", design = 7L, target = "pleiotropy", t = 1L),
  list(scen = "rp", design = 9L, target = "replication", t = 1L),
  list(scen = "med", design = 6L, target = "mediation", t = 1L)
)

for (case in cases) {
  K <- 4L
  set.seed(sim_seed(11L, case$design, K, case$t, 1L, 1L))
  generated <- make_dataset(case$scen, case$design, K, case$t, 1L, N = N)$train
  z <- generated$Z
  estimated <- fit_case(z, K, FALSE)
  stopifnot(!is.null(estimated$params))
  raw <- csm_lfdr(estimated, z, case$target, case$t, cores = cores)
  raw_rejected <- csm_stepup(raw, 0.1)$rejected
  raw_score <- incon_score(z, raw, case$target, case$t)
  raw_decision <- incon_decision(z, raw_rejected, case$target, case$t)

  adjusted <- csm_adjust(z, raw, case$target, case$t, q = 0.1, all = TRUE)
  adjusted_score <- incon_score(z, adjusted$adjusted, case$target, case$t)
  adjusted_decision <- incon_decision(z, adjusted$rejected, case$target, case$t)

  independent <- fit_case(z, K, TRUE)
  independent_lfdr <- csm_lfdr(independent, z, case$target, case$t, cores = cores)
  independent_score <- incon_score(z, independent_lfdr, case$target, case$t)

  t_boundary <- K - 1L
  boundary_lfdr <- csm_lfdr(estimated, z, case$target, t_boundary, cores = cores)
  boundary_score <- incon_score(z, boundary_lfdr, case$target, t_boundary)

  say(sprintf(paste0("%s: raw score=%s, raw decision=%d; majorant score=%s, ",
                     "majorant decision=%d; fixed-Sigma score=%s; t=K-1 score=%s"),
              case$target, raw_score, raw_decision,
              adjusted_score, adjusted_decision,
              independent_score, boundary_score))
  stopifnot(adjusted_score == 0L, adjusted_decision == 0L)

  if (requireNamespace("csmGmm", quietly = TRUE) && K <= 3L) {
    upstream <- csmGmm::check_incongruous(z, raw)
    say("upstream orthant diagnostic count=", length(upstream))
  } else {
    say("upstream orthant diagnostic skipped: it supports only K=2 or K=3")
  }
}

say("DIAGNOSIS_DONE")
