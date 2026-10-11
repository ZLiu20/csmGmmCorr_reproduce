minimax_pval_R1 <- function(t) {
  a <- abs(t)
  M <- pmax(a[, 1L], a[, 2L])
  m <- pmin(a[, 1L], a[, 2L])
  uM <- stats::pnorm(M, lower.tail = FALSE)
  um <- stats::pnorm(m, lower.tail = FALSE)
  gap <- um - uM
  K <- ifelse(gap > 0, floor(uM / gap), Inf)
  K <- pmin(K, 1e12)
  harmonic <- digamma(K + 1) - digamma(1)
  pmin(1, 2 * um / (K + 1) + 2 * harmonic * gap)
}

run_mediation_test_R1 <- function(z, t_value, q_value) {
  if (!requireNamespace("mediation.test", quietly = TRUE)) stop("The mediation.test package is required")
  z <- as.matrix(z)
  if (ncol(z) %% 2L) stop("Mediation requires even K")
  pairs <- matrix(seq_len(ncol(z)), ncol = 2L, byrow = TRUE)
  pair_p <- matrix(1, nrow(z), nrow(pairs))
  for (j in seq_len(nrow(pairs))) {
    mediation.test::mediation_test_minimax(t = z[, pairs[j, ], drop = FALSE], alpha = q_value,
                                           truncation = 0, sample_size = Inf, compute_pvals = FALSE)
    pair_p[, j] <- clamp_probability_R1(minimax_pval_R1(z[, pairs[j, ], drop = FALSE]))
  }
  combined <- combine_mediation_target_pvalues_R1(z, pair_p, t_value)
  adjusted <- stats::p.adjust(combined$global, method = "BH")
  list(rejected = adjusted <= q_value, score = combined$global,
       rankingScoreType = "mediation_intersection_global_p",
       tuning = paste0("truncation=0;sample_size=Inf;pathways=", nrow(pairs),
                       ";support_partial_conjunction_r=", as.integer(t_value) + 1L,
                       ";pathway_any=Bonferroni;combine=max;FDR=BH"))
}
