# MediationTest (minimax test, package mediation.test): mediation.
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
      compute_pvals = TRUE
    )
    pair_p[, j] <- clamp_probability_R1(fit$pvals)
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

