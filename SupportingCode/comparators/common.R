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
    if (ncol(beta) %% 2L) stop("Mediation requires an even number of coordinates")
    pair_mat <- matrix(seq_len(ncol(beta)), ncol = 2L, byrow = TRUE)
  }
  complete <- matrix(FALSE, nrow(beta), nrow(pair_mat))
  for (j in seq_len(nrow(pair_mat))) {
    complete[, j] <- beta[, pair_mat[j, 1L]] != 0 & beta[, pair_mat[j, 2L]] != 0
  }
  t_value <- as.integer(t_value)
  if (!is.finite(t_value) || t_value < 1L || t_value >= ncol(beta)) {
    stop("For mediation, t must satisfy 1 <= t < K")
  }
  support > t_value & rowSums(complete) >= 1L
}

score_comparison_decisions_R1 <- function(rejected, score, beta, analysis,
                                          t_value, pair_mat = NULL) {
  rejected <- as.logical(rejected)
  score <- clamp_probability_R1(score)
  beta <- as.matrix(beta)
  if (length(rejected) != nrow(beta) || length(score) != nrow(beta)) {
    stop("Decision, score, and truth vectors have incompatible lengths")
  }
  alt <- truth_for_comparison_R1(beta, analysis, t_value, pair_mat)
  support <- rowSums(beta != 0)
  structural_null <- if (analysis == "mediation") support > as.integer(t_value) & !alt else rep(FALSE, nrow(beta))
  n_rej <- sum(rejected)
  n_alt <- sum(alt)
  tp <- sum(rejected & alt)
  fp <- sum(rejected & !alt)
  list(
    summary = data.frame(
      nAlternative = n_alt,
      nRej = n_rej,
      truePositive = tp,
      falsePositive = fp,
      structuralNullFalsePositive = sum(rejected & structural_null),
      FDP = if (n_rej) fp / n_rej else 0,
      power = if (n_alt) tp / n_alt else 0,
      stringsAsFactors = FALSE
    ),
    rejected = rejected,
    score = score,
    truthAlternative = alt
  )
}

combine_pathway_pvalues_R1 <- function(pair_pvalues, r = 1L) {
  pair_pvalues <- as.matrix(pair_pvalues)
  pair_pvalues[!is.finite(pair_pvalues)] <- 1
  pair_pvalues[] <- pmin(1, pmax(0, pair_pvalues))
  paths <- ncol(pair_pvalues)
  r <- as.integer(r)
  if (!is.finite(r) || r < 1L || r > paths) stop("r must index an existing pathway")
  ordered <- t(apply(pair_pvalues, 1L, sort))
  pmin(1, (paths - r + 1L) * ordered[, r])
}

coordinate_partial_conjunction_p_R1 <- function(z, t_value) {
  p <- two_sided_p_R1(z)
  K <- ncol(p)
  t_value <- as.integer(t_value)
  if (!is.finite(t_value) || t_value < 1L || t_value >= K) stop("t must satisfy 1 <= t < K")
  ordered <- t(apply(p, 1L, sort))
  pmin(1, (K - t_value) * ordered[, t_value + 1L])
}

combine_any_pathway_pvalues_R1 <- function(pair_pvalues) {
  pair_pvalues <- as.matrix(pair_pvalues)
  pair_pvalues[!is.finite(pair_pvalues)] <- 1
  pair_pvalues[] <- pmin(1, pmax(0, pair_pvalues))
  pmin(1, ncol(pair_pvalues) * apply(pair_pvalues, 1L, min))
}

combine_mediation_target_pvalues_R1 <- function(z, pair_pvalues, t_value) {
  coordinate <- coordinate_partial_conjunction_p_R1(z, t_value)
  pathway <- combine_any_pathway_pvalues_R1(pair_pvalues)
  list(global = pmax(coordinate, pathway), coordinate = coordinate, pathway = pathway)
}
