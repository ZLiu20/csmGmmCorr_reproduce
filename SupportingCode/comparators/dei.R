dei_b_pvalue_R1 <- function(px, py) {
  px <- clamp_probability_R1(px)
  py <- clamp_probability_R1(py)
  if (length(px) != length(py)) stop("DEI inputs differ in length")
  n <- length(px)
  pmax_pair <- pmax(px, py)
  ordering <- order(pmax_pair, seq_along(pmax_pair))
  sorted <- pmax_pair[ordering]
  estimate <- rep(NA_real_, n)
  previous <- 0
  local_min <- 1L
  for (i in seq_len(n)) {
    epsilon <- sorted[i]
    if (!is.finite(epsilon) || epsilon >= 0.4) break
    while (local_min <= n && sorted[local_min] <= epsilon) local_min <- local_min + 1L
    region <- local_min - i
    current <- (1 - 2 * epsilon) * (region / n - epsilon)
    current <- max(previous, current, 0)
    estimate[i] <- min(1, current)
    previous <- estimate[i]
  }
  out <- rep(1, n)
  out[ordering] <- clamp_probability_R1(estimate)
  out
}

run_dei_mediation_R1 <- function(z, t_value, q_value) {
  z <- as.matrix(z)
  if (ncol(z) %% 2L) stop("Mediation requires even K")
  p <- two_sided_p_R1(z)
  pairs <- matrix(seq_len(ncol(z)), ncol = 2L, byrow = TRUE)
  pair_p <- matrix(1, nrow(z), nrow(pairs))
  for (j in seq_len(nrow(pairs))) pair_p[, j] <- dei_b_pvalue_R1(p[, pairs[j, 1L]], p[, pairs[j, 2L]])
  combined <- combine_mediation_target_pvalues_R1(z, pair_p, t_value)
  adjusted <- stats::p.adjust(combined$global, method = "BH")
  list(rejected = adjusted <= q_value, score = combined$global,
       rankingScoreType = "mediation_intersection_global_p",
       tuning = paste0("pathways=", nrow(pairs), ";support_partial_conjunction_r=",
                       as.integer(t_value) + 1L, ";pathway_any=Bonferroni;combine=max;FDR=BH"))
}
