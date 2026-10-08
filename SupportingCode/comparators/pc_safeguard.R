# Parameter-free partial-conjunction safeguard: mediation.
run_mediation_safeguard_R1 <- function(z, t_value, q_value) {
  z <- as.matrix(z)
  if (ncol(z) %% 2L != 0L) stop("Mediation requires even K.")
  p_matrix <- two_sided_p_R1(z)
  pair_mat <- matrix(seq_len(ncol(z)), ncol = 2L, byrow = TRUE)
  pair_p <- matrix(1, nrow(z), nrow(pair_mat))
  for (j in seq_len(nrow(pair_mat))) {
    pair_p[, j] <- pmax(
      p_matrix[, pair_mat[j, 1L]], p_matrix[, pair_mat[j, 2L]]
    )
  }
  combined <- combine_mediation_target_pvalues_R1(z, pair_p, t_value)
  global_p <- combined$global
  adjusted <- stats::p.adjust(global_p, method = "BH")
  list(
    rejected = adjusted <= q_value,
    score = global_p,
    rankingScoreType = "parameter_free_mediation_safeguard_p",
    tuning = paste0(
      "support_partial_conjunction_r=", as.integer(t_value) + 1L,
      ";pair=max_raw_coordinate_p;pathway_any=Bonferroni;combine=max;FDR=BH"
    )
  )
}

