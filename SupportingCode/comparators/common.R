# Shared helpers for the comparison methods (p-values, truth, mediation combination rules).
clamp_probability_R1 <- function(x, replacement = 1) {
  x <- as.numeric(x)
  x[!is.finite(x)] <- replacement
  pmin(1, pmax(0, x))
}

two_sided_p_R1 <- function(z) {
  p <- 2 * stats::pnorm(-abs(as.matrix(z)))
  p[p < .Machine$double.xmin] <- .Machine$double.xmin
  p[p > 1 - .Machine$double.eps] <- 1 - .Machine$double.eps
  p
}

truth_for_comparison_R1 <- function(beta, analysis, t_value, pair_mat = NULL) {
  beta <- as.matrix(beta)
  support <- rowSums(beta != 0)
  if (analysis == "pleiotropy") return(support > t_value)
  if (analysis == "replication") {
    same_direction <- support == abs(rowSums(sign(beta)))
    return(support > t_value & same_direction)
  }
  if (analysis != "mediation") stop("Unknown comparison analysis: ", analysis)
  if (is.null(pair_mat)) {
    if (ncol(beta) %% 2L != 0L) stop("Mediation requires an even K.")
    pair_mat <- matrix(seq_len(ncol(beta)), ncol = 2L, byrow = TRUE)
  }
  complete <- matrix(FALSE, nrow(beta), nrow(pair_mat))
  for (j in seq_len(nrow(pair_mat))) {
    complete[, j] <- beta[, pair_mat[j, 1L]] != 0 & beta[, pair_mat[j, 2L]] != 0
  }
  t_value <- as.integer(t_value)
  if (!is.finite(t_value) || t_value < 1L || t_value >= ncol(beta)) {
    stop("For mediation, t_value must satisfy 1 <= t < K.")
  }
  support > t_value & rowSums(complete) >= 1L
}

score_comparison_decisions_R1 <- function(rejected, score, beta, analysis,
                                          t_value, pair_mat = NULL) {
  rejected <- as.logical(rejected)
  score <- clamp_probability_R1(score)
  beta <- as.matrix(beta)
  if (length(rejected) != nrow(beta) || length(score) != nrow(beta)) {
    stop("Comparison decisions, score, and beta have incompatible lengths.")
  }
  truth_alt <- truth_for_comparison_R1(beta, analysis, t_value, pair_mat)
  support <- rowSums(beta != 0)
  structural_null <- if (analysis == "mediation") {
    support > as.integer(t_value) & !truth_alt
  } else {
    rep(FALSE, nrow(beta))
  }
  n_rej <- sum(rejected)
  n_alt <- sum(truth_alt)
  tp <- sum(rejected & truth_alt)
  fp <- sum(rejected & !truth_alt)
  list(
    summary = data.frame(
      nAlternative = n_alt,
      nRej = n_rej,
      truePositive = tp,
      falsePositive = fp,
      structuralNullFalsePositive = sum(rejected & structural_null),
      FDP = if (n_rej > 0L) fp / n_rej else 0,
      Power = if (n_alt > 0L) tp / n_alt else 0,
      stringsAsFactors = FALSE
    ),
    rejected = rejected,
    score = score,
    truthAlternative = truth_alt
  )
}

combine_pathway_pvalues_R1 <- function(pair_pvalues, r = 1L) {
  pair_pvalues <- as.matrix(pair_pvalues)
  pair_pvalues[!is.finite(pair_pvalues)] <- 1
  pair_pvalues[pair_pvalues < 0] <- 0
  pair_pvalues[pair_pvalues > 1] <- 1
  n_paths <- ncol(pair_pvalues)
  r <- as.integer(r)
  if (!is.finite(r) || r < 1L || r > n_paths) {
    stop("The pathway partial-conjunction order r must lie in 1..number of pathways.")
  }
  ordered <- t(apply(pair_pvalues, 1L, sort))
  # Bonferroni partial conjunction: test that at least r pathway nulls are false.
  pmin(1, (n_paths - r + 1L) * ordered[, r])
}

coordinate_partial_conjunction_p_R1 <- function(z, t_value) {
  p_matrix <- two_sided_p_R1(z)
  K <- ncol(p_matrix)
  t_value <- as.integer(t_value)
  if (!is.finite(t_value) || t_value < 1L || t_value >= K) {
    stop("The coordinate threshold t_value must satisfy 1 <= t < K.")
  }
  ordered <- t(apply(p_matrix, 1L, sort))
  pmin(1, (K - t_value) * ordered[, t_value + 1L])
}

combine_any_pathway_pvalues_R1 <- function(pair_pvalues) {
  pair_pvalues <- as.matrix(pair_pvalues)
  pair_pvalues[!is.finite(pair_pvalues)] <- 1
  pair_pvalues[pair_pvalues < 0] <- 0
  pair_pvalues[pair_pvalues > 1] <- 1
  pmin(1, ncol(pair_pvalues) * apply(pair_pvalues, 1L, min))
}

# Mediation null at threshold t: at most t nonzero links, or no complete pathway. Its p-value is
# the larger of the coordinate partial-conjunction p-value and the pathway Bonferroni p-value.
# At t = 1 a complete pathway already has two nonzero links, so the pathway p-value alone is
# valid for the whole null and the partial-conjunction term is omitted (it can only enlarge the
# p-value). options(csm.med_pc_screen_t1 = TRUE) restores the earlier max() at t = 1 (testing).
med_combine_rule_R1 <- function(t_value)
  if (as.integer(t_value) == 1L && !isTRUE(getOption("csm.med_pc_screen_t1", FALSE))) "pathway" else "max"
combine_mediation_target_pvalues_R1 <- function(z, pair_pvalues, t_value) {
  p_coordinate <- coordinate_partial_conjunction_p_R1(z, t_value)
  p_pathway <- combine_any_pathway_pvalues_R1(pair_pvalues)
  list(
    global = if (med_combine_rule_R1(t_value) == "pathway") p_pathway else pmax(p_coordinate, p_pathway),
    coordinate = p_coordinate,
    pathway = p_pathway
  )
}

