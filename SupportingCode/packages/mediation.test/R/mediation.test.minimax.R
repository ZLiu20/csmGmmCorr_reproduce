mediation_test_minimax <- function(t, alpha = 0.05, truncation = 0,
                                   sample_size = Inf, compute_pvals = TRUE) {
  t <- .mt_matrix(t)
  alpha <- .mt_scalar(alpha, 0, 0.5, "alpha")
  truncation <- .mt_scalar(truncation, 0, Inf, "truncation")
  sample_size <- .mt_sample_size(sample_size)
  compute_pvals <- .mt_flag(compute_pvals, "compute_pvals")

  decide <- function(level) {
    inverse_level <- floor(1 / level)
    grid <- seq(
      1 - inverse_level * level / 2,
      1 - inverse_level * level / 2 + inverse_level * level / 2,
      by = level / 2
    )
    cutoffs <- .mt_quantile(grid, sample_size)
    absolute_t <- abs(t)
    answer <- findInterval(absolute_t[, 1L], cutoffs) ==
      findInterval(absolute_t[, 2L], cutoffs)
    names(answer) <- rep("rejection", length(answer))
    answer
  }

  conservative_pvalues <- function() {
    absolute_t <- abs(t)
    upper_probability <- .mt_cdf(
      pmax(absolute_t[, 1L], absolute_t[, 2L]), sample_size
    )
    lower_probability <- .mt_cdf(
      pmin(absolute_t[, 1L], absolute_t[, 2L]), sample_size
    )
    gap <- upper_probability - lower_probability

    number_of_steps <- rep(Inf, length(gap))
    separated <- gap > 0
    number_of_steps[separated] <- floor(
      (1 - upper_probability[separated]) / gap[separated]
    )

    harmonic <- numeric(length(gap))
    finite_steps <- is.finite(number_of_steps) & number_of_steps > 0
    harmonic[finite_steps] <- digamma(number_of_steps[finite_steps] + 1) +
      0.5772156649015329

    pvalue <- numeric(length(gap))
    pvalue[separated] <-
      2 * (1 - lower_probability[separated]) /
        (number_of_steps[separated] + 1) +
      2 * gap[separated] * harmonic[separated]
    pmin(1, pmax(0, pvalue))
  }

  decision <- decide(alpha)
  if (truncation > 0) {
    decision <- decision & apply(abs(t), 1L, min) >= truncation
    names(decision) <- rep("rejection", length(decision))
  }

  pval <- if (compute_pvals) conservative_pvalues() else NA_real_
  answer <- list(
    t = t,
    alpha = alpha,
    truncation = truncation,
    sample_size = sample_size,
    decision = decision,
    pval = pval,
    method = "minimax"
  )
  class(answer) <- "mediation.test"
  answer
}

mediation_test_minimax_BH <- function(t, alpha = 0.05, truncation = 0,
                                      BH = c("statistics", "pval"),
                                      sample_size = Inf, K = 10L) {
  t <- .mt_matrix(t)
  alpha <- .mt_scalar(alpha, 0, 0.5, "alpha")
  truncation <- .mt_scalar(truncation, 0, Inf, "truncation")
  sample_size <- .mt_sample_size(sample_size)
  K <- .mt_integer(K, 1L, Inf, "K")
  BH <- match.arg(BH)
  number_of_tests <- nrow(t)

  if (number_of_tests == 1L) {
    warning("One single test statistic. Switching to 'mediation_test_minimax'.")
    return(mediation_test_minimax(t, alpha, truncation, sample_size))
  }
  if (K >= number_of_tests) {
    K <- number_of_tests - 1L
    warning("Replacing the user-supplied 'K' with 'nrow(t)-1'.")
  }

  if (BH == "pval") {
    answer <- mediation_test_minimax(
      t, alpha, truncation, sample_size, compute_pvals = TRUE
    )
    adjusted_pvalue <- stats::p.adjust(answer$pval)
    answer$decision <- adjusted_pvalue <= alpha &
      apply(abs(t), 1L, min) >= truncation
    names(answer$decision) <- rep("rejection", number_of_tests)
  } else {
    number_rejected <- function(level) {
      sum(mediation_test_minimax(
        t, level, truncation, sample_size, compute_pvals = FALSE
      )$decision)
    }

    keep_searching <- TRUE
    from <- -K + 1L
    to <- 0L
    chosen_level <- NA_real_
    valid_levels <- logical(0L)
    cumulative_valid <- integer(0L)

    while (keep_searching && to < number_of_tests) {
      from <- from + K
      to <- min(to + K, number_of_tests)
      indices <- seq.int(from, to)
      candidate_levels <- alpha * indices / number_of_tests
      rejected <- vapply(candidate_levels, number_rejected, integer(1L))
      valid <- rejected >= indices
      valid_levels <- c(valid_levels, valid)
      cumulative_valid <- cumsum(valid_levels)

      earlier <- if (length(cumulative_valid) > K) {
        cumulative_valid[seq_len(length(cumulative_valid) - K)]
      } else {
        integer(0L)
      }
      lagged <- c(rep(NA_integer_, K), earlier)
      stalled <- any(cumulative_valid - lagged <= 0, na.rm = TRUE)
      if (stalled) {
        keep_searching <- FALSE
        if (max(cumulative_valid, na.rm = TRUE) == 0L) {
          chosen_level <- alpha
          answer <- mediation_test_minimax(
            t, chosen_level, truncation, sample_size,
            compute_pvals = FALSE
          )
          answer$decision[] <- FALSE
        } else {
          final_index <- which.max(cumulative_valid)
          chosen_level <- alpha * final_index / number_of_tests
          answer <- mediation_test_minimax(
            t, chosen_level, truncation, sample_size,
            compute_pvals = FALSE
          )
          answer$alpha <- alpha
        }
      }
    }

    if (is.na(chosen_level)) {
      answer <- mediation_test_minimax(
        t, alpha, truncation, sample_size, compute_pvals = FALSE
      )
    }
  }

  answer$method <- "minimax_BH"
  answer$BH <- BH
  answer$K <- K
  answer
}
