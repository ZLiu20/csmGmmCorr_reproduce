adafilter_bh_R1 <- function(p_matrix, r, alpha = 0.1) {
  p_matrix <- as.matrix(p_matrix)
  p_matrix[!is.finite(p_matrix)] <- NA_real_
  available <- rowSums(!is.na(p_matrix))
  eligible <- which(available >= r)
  answer <- data.frame(
    decision = rep(FALSE, nrow(p_matrix)), adjusted.p = NA_real_,
    selection.p = NA_real_, filter.p = NA_real_, adj.number = NA_real_
  )
  if (!length(eligible)) return(answer)
  ordered <- t(apply(p_matrix[eligible, , drop = FALSE], 1L, sort, na.last = TRUE))
  selection <- ordered[, r] * (available[eligible] - r + 1L)
  filter <- ordered[, 1L] * (available[eligible] - r + 1L)
  ordering <- order(selection, seq_along(selection))
  sorted <- selection[ordering]
  adjustment_number <- vapply(seq_along(sorted), function(i) sum(filter <= sorted[i], na.rm = TRUE), numeric(1))
  adjusted_sorted <- adjustment_number / length(sorted) * sorted
  adjusted_sorted <- pmin(1, rev(cummin(rev(adjusted_sorted))))
  adjusted <- numeric(length(sorted))
  adjusted[ordering] <- adjusted_sorted
  original_number <- numeric(length(sorted))
  original_number[ordering] <- adjustment_number
  answer$decision[eligible] <- adjusted <= alpha
  answer$adjusted.p[eligible] <- adjusted
  answer$selection.p[eligible] <- selection
  answer$filter.p[eligible] <- filter
  answer$adj.number[eligible] <- original_number
  answer
}

run_adafilter_pleiotropy_R1 <- function(z, t_value, q_value) {
  p <- two_sided_p_R1(z)
  r <- as.integer(t_value) + 1L
  result <- if (requireNamespace("adaFilter", quietly = TRUE)) {
    adaFilter::adaFilter(p.matrix = p, r = r, type.I.err = "FDR", alpha = q_value, fast = FALSE)
  } else adafilter_bh_R1(p, r, q_value)
  if (!is.data.frame(result) || nrow(result) != nrow(p)) stop("Invalid adaFilter result")
  adjusted <- clamp_probability_R1(result$adjusted.p)
  list(rejected = as.logical(result$decision %in% TRUE), score = adjusted,
       rankingScoreType = "adjusted_p", tuning = paste0("r=", r, ";type=FDR"))
}
