compute_map_rejection_probs <- function(alpha = 0.05, K = 16L,
                                        loss = c("0-1", "quadratic"),
                                        sample_size = Inf, truncation = 0,
                                        return_solver = FALSE) {
  alpha <- .mt_scalar(alpha, 0, 0.5, "alpha")
  K <- .mt_integer(K, 2L, 128L, "K")
  loss <- match.arg(loss)
  sample_size <- .mt_sample_size(sample_size)
  truncation <- .mt_scalar(truncation, 0, Inf, "truncation")
  return_solver <- .mt_flag(return_solver, "return_solver")

  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("Package 'lpSolve' is required to compute a rejection map.", call. = FALSE)
  }

  boundary <- 2 * .mt_quantile(1 - alpha / 2, sample_size)
  prior_sd <- sqrt(boundary)
  cell_edges <- seq(-boundary, boundary, length.out = 2L * K + 1L)
  null_grid <- seq(-2 * boundary, 2 * boundary, length.out = 4L * K + 1L)
  null_x <- c(null_grid, rep(0, 4L * K))
  null_y <- c(rep(0, 4L * K + 1L), null_grid[-(2L * K + 1L)])

  outside_probability <- function(location, cutoff) {
    1 - stats::pt(cutoff - location, df = sample_size - 1L) +
      stats::pt(-cutoff - location, df = sample_size - 1L)
  }
  upper_tail <- outside_probability(null_x, boundary / 2) *
    outside_probability(null_y, boundary / 2)
  corner_adjustment <-
    (outside_probability(null_x, boundary / 2) -
       outside_probability(null_x, boundary)) *
    (outside_probability(null_y, boundary / 2) -
       outside_probability(null_y, boundary))
  constraint_rhs <- alpha - upper_tail + corner_adjustment

  cell_allowed <- rep(1, 4L * K^2)
  if (truncation > 0) {
    away_from_zero <- cell_edges[-length(cell_edges)] > truncation |
      cell_edges[-1L] < -truncation
    cell_allowed <- as.vector(outer(away_from_zero, away_from_zero, "*"))
  }
  constraint_rhs <- c(constraint_rhs, cell_allowed)
  constraint_direction <- rep("<=", length(constraint_rhs))

  interval_weights <- array(0, dim = c(2L * K, 2L * K, length(null_x)))
  for (i in seq_len(2L * K)) {
    probability_x <- stats::pt(
      cell_edges[i + 1L] - null_x, df = sample_size - 1L
    ) - stats::pt(cell_edges[i] - null_x, df = sample_size - 1L)
    for (j in seq_len(2L * K)) {
      probability_y <- stats::pt(
        cell_edges[j + 1L] - null_y, df = sample_size - 1L
      ) - stats::pt(cell_edges[j] - null_y, df = sample_size - 1L)
      interval_weights[i, j, ] <- probability_x * probability_y
    }
  }
  constraint_matrix <- cbind(
    apply(interval_weights, 3L, as.vector),
    diag(4L * K^2)
  )

  marginal_risk <- function(power = 0L) {
    if (is.infinite(sample_size)) {
      total_sd <- sqrt(1 + prior_sd^2)
      upper <- cell_edges[-1L]
      lower <- cell_edges[-length(cell_edges)]
      if (power == 0L) {
        return(stats::pnorm(upper / total_sd) - stats::pnorm(lower / total_sd))
      }
      first_term <- prior_sd^2 *
        (stats::pnorm(upper / total_sd) - stats::pnorm(lower / total_sd))
      second_term <- prior_sd^4 / (1 + prior_sd^2)^(3 / 2) *
        (upper * stats::dnorm(upper / total_sd) -
           lower * stats::dnorm(lower / total_sd))
      return(first_term - second_term)
    }

    integrated <- vapply(cell_edges, function(upper) {
      integrand <- function(location) {
        location^power * stats::pt(
          upper, ncp = location, df = sample_size - 1L
        ) * stats::dnorm(location, sd = prior_sd)
      }
      stats::integrate(integrand, -Inf, Inf)$value
    }, numeric(1L))
    diff(integrated)
  }

  one_dimensional_risk <- marginal_risk(if (loss == "0-1") 0L else 2L)
  risk <- outer(one_dimensional_risk, one_dimensional_risk, "*")
  solver <- lpSolve::lp(
    direction = "max",
    objective.in = as.vector(risk),
    const.mat = constraint_matrix,
    const.dir = constraint_direction,
    const.rhs = constraint_rhs,
    transpose.constraints = FALSE
  )
  if (!identical(solver$status, 0L)) {
    stop(sprintf("The map linear program failed with status %s.", solver$status),
         call. = FALSE)
  }

  map <- matrix(
    pmin(1, pmax(0, solver$solution)),
    nrow = 2L * K,
    ncol = 2L * K
  )
  attr(map, "solver") <- if (return_solver) solver else NA
  attr(map, "loss") <- loss
  attr(map, "alpha") <- alpha
  attr(map, "truncation") <- truncation
  attr(map, "sample_size") <- sample_size
  class(map) <- "map.rejection.probabilities"
  map
}
