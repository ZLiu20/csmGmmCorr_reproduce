# adaFilter (Wang et al.): pleiotropy, partial conjunction r = t + 1.
# Adaptive-filtering BH for an r-of-K partial conjunction hypothesis. This is
# the algorithm used by adaFilter 1.0.0 and keeps the code runnable when the
# optional package is not installed on a compute node.
adafilter_bh_R1 <- function(p_matrix, r, alpha = 0.1) {
  p_matrix <- as.matrix(p_matrix)
  p_matrix[!is.finite(p_matrix)] <- NA_real_
  n_hypotheses <- nrow(p_matrix)
  n_available <- rowSums(!is.na(p_matrix))
  eligible <- which(n_available >= r)
  result <- data.frame(
    decision = rep(FALSE, n_hypotheses),
    adjusted.p = rep(NA_real_, n_hypotheses),
    selection.p = rep(NA_real_, n_hypotheses),
    filter.p = rep(NA_real_, n_hypotheses),
    adj.number = rep(NA_real_, n_hypotheses)
  )
  if (length(eligible) == 0L) return(result)

  p_sorted <- t(apply(p_matrix[eligible, , drop = FALSE], 1L, sort, na.last = TRUE))
  n_available_sub <- n_available[eligible]
  selection <- p_sorted[, r] * (n_available_sub - r + 1L)
  filter <- p_sorted[, r - 1L] * (n_available_sub - r + 1L)
  n <- length(eligible)
  selection_order <- order(selection, seq_along(selection))
  selection_sorted <- selection[selection_order]

  # A stable rank construction avoids relying on duplicate vector names.
  adjustment_number <- vapply(seq_len(n), function(i) {
    sum(filter <= selection_sorted[i], na.rm = TRUE)
  }, numeric(1L))
  adjusted_sorted <- adjustment_number / seq_len(n) * selection_sorted
  adjusted_sorted <- pmin(1, rev(cummin(rev(adjusted_sorted))))
  adjusted <- numeric(n)
  adjustment_original <- numeric(n)
  adjusted[selection_order] <- adjusted_sorted
  adjustment_original[selection_order] <- adjustment_number

  result$decision[eligible] <- adjusted <= alpha
  result$adjusted.p[eligible] <- adjusted
  result$selection.p[eligible] <- selection
  result$filter.p[eligible] <- filter
  result$adj.number[eligible] <- adjustment_original
  result
}

run_adafilter_pleiotropy_R1 <- function(z, t_value, q_value) {
  p_matrix <- two_sided_p_R1(z)
  r <- as.integer(t_value) + 1L
  result <- if (requireNamespace("adaFilter", quietly = TRUE)) {
    adaFilter::adaFilter(
      p.matrix = p_matrix, r = r, type.I.err = "FDR",
      alpha = q_value, fast = FALSE
    )
  } else {
    adafilter_bh_R1(p_matrix, r = r, alpha = q_value)
  }
  if (!is.data.frame(result) || nrow(result) != nrow(p_matrix)) {
    stop("adaFilter returned an invalid result.")
  }
  adjusted <- clamp_probability_R1(result$adjusted.p)
  list(
    rejected = as.logical(result$decision %in% TRUE),
    score = adjusted,
    rankingScoreType = "adjusted_p",
    tuning = paste0("r=", r, ";type=FDR")
  )
}

