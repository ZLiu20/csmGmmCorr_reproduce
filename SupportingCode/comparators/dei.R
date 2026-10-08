# DEI (Liu, Biometrics 2025): mediation, pathway p-values combined over mediators.
# GPLv3 DEI_B algorithm (Liu, Biometrics 2025), implemented here so the method
# job does not depend on a non-CRAN package installation on every node.
dei_b_pvalue_R1 <- function(px, py) {
  px <- clamp_probability_R1(px)
  py <- clamp_probability_R1(py)
  n <- length(px)
  if (length(py) != n) stop("DEI input p-value vectors differ in length.")
  p_min_sorted <- sort(pmin(px, py))
  p_max <- pmax(px, py)
  max_order <- order(p_max, seq_along(p_max))
  p_max_sorted <- p_max[max_order]
  estimated_sorted <- rep(NA_real_, n)
  previous <- 0
  loc_min <- 1L
  for (i in seq_len(n)) {
    epsilon <- p_max_sorted[i]
    if (!is.finite(epsilon) || epsilon >= 0.4) break
    while (loc_min <= n && p_min_sorted[loc_min] <= epsilon) {
      loc_min <- loc_min + 1L
    }
    n_region <- loc_min - i
    estimate <- epsilon / (1 - 2 * epsilon) * (n_region / n - epsilon)
    estimate <- max(previous, estimate, 0)
    estimated_sorted[i] <- min(1, estimate)
    previous <- estimated_sorted[i]
  }
  out <- rep(1, n)
  out[max_order] <- clamp_probability_R1(estimated_sorted)
  out
}

run_dei_mediation_R1 <- function(z, t_value, q_value) {
  z <- as.matrix(z)
  if (ncol(z) %% 2L != 0L) stop("Mediation requires even K.")
  p_matrix <- two_sided_p_R1(z)
  pair_mat <- matrix(seq_len(ncol(z)), ncol = 2L, byrow = TRUE)
  pair_p <- matrix(1, nrow(z), nrow(pair_mat))
  for (j in seq_len(nrow(pair_mat))) {
    pair_p[, j] <- dei_b_pvalue_R1(
      p_matrix[, pair_mat[j, 1L]], p_matrix[, pair_mat[j, 2L]]
    )
  }
  combined <- combine_mediation_target_pvalues_R1(z, pair_p, t_value)
  global_p <- combined$global
  adjusted <- stats::p.adjust(global_p, method = "BH")
  list(
    rejected = adjusted <= q_value,
    score = global_p,
    rankingScoreType = "mediation_intersection_global_p",
    tuning = paste0(
      "pathways=", nrow(pair_mat),
      ";support_partial_conjunction_r=", as.integer(t_value) + 1L,
      ";pathway_any=Bonferroni;combine=max;FDR=BH"
    )
  )
}

