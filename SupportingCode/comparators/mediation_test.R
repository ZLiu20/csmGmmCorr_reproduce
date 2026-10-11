# MediationTest (minimax test, package mediation.test): mediation.
# Conservative p-value of the minimax test (Miles and Chambaz), as in mediation.test with
# sample_size = Inf: p = 2 (1 - Phi(m)) / (K + 1) + 2 (Phi(M) - Phi(m)) H_K, m = min |t|, M = max |t|,
# K = floor((1 - Phi(M)) / (Phi(M) - Phi(m))), H_K the K-th harmonic number. The package builds
# the vector 1:max(K), which fails (NA/NaN) when |t1| = |t2| (e.g. statistics capped at the same
# value) and is very large when they are nearly equal; here H_K = digamma(K + 1) + Euler's constant
# and upper tails are computed directly, which gives the same values otherwise.
minimax_pval_R1 <- function(t) {
  a <- abs(t); M <- pmax(a[, 1], a[, 2]); m <- pmin(a[, 1], a[, 2])
  uM <- stats::pnorm(M, lower.tail = FALSE); um <- stats::pnorm(m, lower.tail = FALSE)
  d <- um - uM
  K <- ifelse(d > 0, floor(uM / d), Inf)
  K <- pmin(K, 1e12)
  H <- digamma(K + 1) - digamma(1)
  pmin(1, 2 * um / (K + 1) + 2 * d * H)
}
run_mediation_test_R1 <- function(z, t_value, q_value) {
  if (!requireNamespace("mediation.test", quietly = TRUE)) {
    stop("The mediation.test package is required for the MediationTest comparator.")
  }
  z <- as.matrix(z)
  if (ncol(z) %% 2L != 0L) stop("Mediation requires even K.")
  pair_mat <- matrix(seq_len(ncol(z)), ncol = 2L, byrow = TRUE)
  pair_p <- matrix(1, nrow(z), nrow(pair_mat))
  for (j in seq_len(nrow(pair_mat))) {
    fit <- mediation.test::mediation_test_minimax(
      t = z[, pair_mat[j, ], drop = FALSE],
      alpha = q_value,
      truncation = 0,
      sample_size = Inf,
      compute_pvals = FALSE
    )
    pair_p[, j] <- clamp_probability_R1(minimax_pval_R1(z[, pair_mat[j, ], drop = FALSE]))
  }
  combined <- combine_mediation_target_pvalues_R1(z, pair_p, t_value)
  global_p <- combined$global
  adjusted <- stats::p.adjust(global_p, method = "BH")
  list(
    rejected = adjusted <= q_value,
    score = global_p,
    rankingScoreType = "mediation_intersection_global_p",
    tuning = paste0(
      "truncation=0;sample_size=Inf;pathways=", nrow(pair_mat),
      ";support_partial_conjunction_r=", as.integer(t_value) + 1L,
      ";pathway_any=Bonferroni;combine=max;FDR=BH"
    )
  )
}

