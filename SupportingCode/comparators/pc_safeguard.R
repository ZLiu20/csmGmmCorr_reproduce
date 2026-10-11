run_mediation_safeguard_R1 <- function(z, t_value, q_value) {
  z <- as.matrix(z)
  if (ncol(z) %% 2L) stop("Mediation requires even K")
  p <- two_sided_p_R1(z)
  pairs <- matrix(seq_len(ncol(z)), ncol = 2L, byrow = TRUE)
  pair_p <- matrix(1, nrow(z), nrow(pairs))
  for (j in seq_len(nrow(pairs))) pair_p[, j] <- pmax(p[, pairs[j, 1L]], p[, pairs[j, 2L]])
  combined <- combine_mediation_target_pvalues_R1(z, pair_p, t_value)
  adjusted <- stats::p.adjust(combined$global, method = "BH")
  list(rejected = adjusted <= q_value, score = combined$global,
       rankingScoreType = "parameter_free_mediation_safeguard_p",
       tuning = paste0("support_partial_conjunction_r=", as.integer(t_value) + 1L,
                       ";pair=max_raw_coordinate_p;pathway_any=Bonferroni;combine=max;FDR=BH"))
}
