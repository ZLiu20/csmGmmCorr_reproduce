mediation_test_Bayes <- function(t, map, sample_size = Inf, truncation = 0) {
  if (!inherits(map, "map.rejection.probabilities")) {
    stop("'map' must be produced by compute_map_rejection_probs().", call. = FALSE)
  }

  truncation <- .mt_scalar(truncation, 0, Inf, "truncation")
  map_truncation <- attr(map, "truncation")
  if (!is.null(map_truncation) && truncation != map_truncation) {
    warning(sprintf(
      "The supplied truncation (%g) differs from the map truncation (%g); using the map value.",
      truncation, map_truncation
    ))
    truncation <- map_truncation
  }

  sample_size <- .mt_sample_size(sample_size)
  map_sample_size <- attr(map, "sample_size")
  if (!is.null(map_sample_size) && sample_size != map_sample_size) {
    warning(sprintf(
      "The supplied sample_size (%s) differs from the map sample_size (%s); using the map value.",
      format(sample_size), format(map_sample_size)
    ))
    sample_size <- map_sample_size
  }

  t <- .mt_matrix(t)
  alpha <- attr(map, "alpha")
  if (is.null(alpha)) stop("The map has no 'alpha' attribute.", call. = FALSE)

  half_grid <- ncol(map) / 2
  boundary <- 2 * .mt_quantile(1 - alpha / 2, sample_size)
  absolute_t <- abs(t)
  largest <- pmax(absolute_t[, 1L], absolute_t[, 2L])
  smallest <- pmin(absolute_t[, 1L], absolute_t[, 2L])

  decision <- largest > boundary & smallest > boundary / 2
  inside_map <- which(largest <= boundary)
  if (length(inside_map)) {
    map_index <- cbind(
      ceiling((t[inside_map, 1L] + boundary) * half_grid / boundary),
      ceiling((t[inside_map, 2L] + boundary) * half_grid / boundary)
    )
    map_index[, 1L] <- pmax(1L, pmin(nrow(map), map_index[, 1L]))
    map_index[, 2L] <- pmax(1L, pmin(ncol(map), map_index[, 2L]))
    decision[inside_map] <- as.logical(map[map_index])
  }
  names(decision) <- rep("rejection", length(decision))

  if (truncation > 0) {
    decision <- decision & apply(absolute_t, 1L, min) >= truncation
    names(decision) <- rep("rejection", length(decision))
  }

  answer <- list(
    t = t,
    alpha = alpha,
    truncation = truncation,
    decision = decision,
    pvals = NA_real_,
    method = "Bayes",
    map = map
  )
  class(answer) <- "mediation.test"
  answer
}
