.mt_matrix <- function(t) {
  if (is.vector(t)) t <- matrix(t, ncol = 2L)
  t <- as.matrix(t)
  storage.mode(t) <- "double"
  if (ncol(t) != 2L || any(!is.finite(t))) {
    stop("'t' must contain finite pairs of test statistics.", call. = FALSE)
  }
  t
}

.mt_scalar <- function(x, lower, upper, name) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      x < lower || x > upper) {
    stop(sprintf("'%s' must be a scalar in [%s, %s].", name, lower, upper),
         call. = FALSE)
  }
  as.numeric(x)
}

.mt_integer <- function(x, lower, upper, name) {
  value <- .mt_scalar(x, lower, upper, name)
  if (!isTRUE(all.equal(value, round(value)))) {
    message(sprintf("Coercing argument '%s' to an integer.", name))
  }
  as.integer(value)
}

.mt_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(sprintf("'%s' must be TRUE or FALSE.", name), call. = FALSE)
  }
  x
}

.mt_quantile <- function(p, sample_size) {
  stats::qt(p, df = sample_size - 1L)
}

.mt_cdf <- function(x, sample_size) {
  stats::pt(x, df = sample_size - 1L)
}

.mt_sample_size <- function(sample_size) {
  value <- .mt_scalar(sample_size, 2, Inf, "sample_size")
  if (is.finite(value) && !isTRUE(all.equal(value, round(value)))) {
    message("Coercing argument 'sample_size' to an integer.")
    value <- as.integer(value)
  }
  value
}

print.mediation.test <- function(x, n = 6L, ...) {
  n <- .mt_integer(n, 1L, Inf, "n")
  shown <- min(n, nrow(x$t))
  method <- x$method

  cat("Testing the composite product-null hypothesis against nonzero products.\n")
  cat("* method:\n")
  if (method == "minimax") {
    cat("minimax\n")
  } else if (method == "minimax_BH") {
    basis <- if (x$BH == "statistics") "test statistics" else "conservative p-values"
    cat(sprintf("minimax multiple testing based on %s\n", basis))
  } else if (method == "Bayes") {
    cat("Bayes-risk optimal\n")
  } else {
    stop("Unknown mediation test method.", call. = FALSE)
  }

  cat("* test statistic:\n")
  print(utils::head(x$t, shown))
  if (nrow(x$t) > shown) cat("...\n")

  if (method == "minimax") {
    cat(sprintf("* wished type-I error: %.3f\n", x$alpha))
  } else if (method == "minimax_BH") {
    cat(sprintf("* wished false discovery rate: %.3f\n", x$alpha))
  } else {
    cat(sprintf("* wished type-I error: %.3f\n", attr(x$map, "alpha")))
    cat(sprintf("* loss function: %s\n", attr(x$map, "loss")))
  }

  cat(sprintf("* user-supplied truncation parameter: %.3f\n", x$truncation))
  if (method == "minimax") {
    cat(sprintf("* sample size used for calibration: %s\n", format(x$sample_size)))
  }
  cat(sprintf("* decision (%d rejection(s) overall):\n", sum(x$decision)))
  text <- ifelse(
    x$decision,
    sprintf("reject at confidence %.3f", 1 - x$alpha),
    sprintf("cannot reject at confidence %.3f", 1 - x$alpha)
  )
  print(utils::head(text, shown))
  if (length(text) > shown) cat("...\n")

  cat("* conservative p-value:\n")
  print(utils::head(x$pval, shown))
  if (length(x$pval) > shown) cat("...\n")
  invisible(x)
}

.mt_require_plotting <- function() {
  needed <- c("ggplot2", "tibble")
  absent <- needed[!vapply(needed, requireNamespace, logical(1L), quietly = TRUE)]
  if (length(absent)) {
    stop(sprintf("Plotting requires: %s.", paste(absent, collapse = ", ")),
         call. = FALSE)
  }
}

plot.mediation.test <- function(x, filename = NULL, return_fig = FALSE,
                                xlim = c(-4, 4), ylim = c(-4, 4), ...) {
  .mt_require_plotting()
  return_fig <- .mt_flag(return_fig, "return_fig")
  if (length(xlim) != 2L || length(ylim) != 2L) {
    stop("'xlim' and 'ylim' must each contain two limits.", call. = FALSE)
  }
  xlim <- sort(as.numeric(xlim))
  ylim <- sort(as.numeric(ylim))
  points <- tibble::tibble(
    statistic_1 = x$t[, 1L],
    statistic_2 = x$t[, 2L],
    accepted = !x$decision
  )

  if (x$method == "Bayes") {
    figure <- plot(x$map, return_fig = TRUE) +
      ggplot2::geom_point(
        data = points,
        ggplot2::aes(x = statistic_1, y = statistic_2, colour = accepted)
      ) +
      ggplot2::ggtitle("Bayes-risk optimal test")
  } else {
    figure <- ggplot2::ggplot(points) +
      ggplot2::geom_point(
        ggplot2::aes(x = statistic_1, y = statistic_2, colour = accepted)
      ) +
      ggplot2::geom_hline(yintercept = 0) +
      ggplot2::geom_vline(xintercept = 0) +
      ggplot2::coord_fixed() +
      ggplot2::labs(
        x = expression(T[1]), y = expression(T[2]),
        title = "Minimax optimal test"
      )
  }
  figure <- figure +
    ggplot2::scale_x_continuous(limits = xlim, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = ylim, expand = c(0, 0)) +
    ggplot2::theme(legend.position = "none")

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, plot = figure)
  } else {
    suppressWarnings(print(figure))
  }
  if (return_fig) figure else invisible(NULL)
}

plot.map.rejection.probabilities <- function(x, filename = NULL,
                                             return_fig = FALSE, ...) {
  .mt_require_plotting()
  return_fig <- .mt_flag(return_fig, "return_fig")
  half_grid <- nrow(x) / 2
  alpha <- attr(x, "alpha")
  sample_size <- attr(x, "sample_size")
  boundary <- 2 * .mt_quantile(1 - alpha / 2, sample_size)
  axis <- seq(-boundary, boundary, length.out = 2L * half_grid)
  grid <- expand.grid(statistic_1 = axis, statistic_2 = axis)
  grid$probability <- as.vector(x)

  subtitle <- if (attr(x, "truncation") > 0) {
    sprintf("type-I error %.3f; truncation %.3f", alpha, attr(x, "truncation"))
  } else {
    sprintf("type-I error %.3f; no truncation", alpha)
  }
  figure <- ggplot2::ggplot(
    grid,
    ggplot2::aes(x = statistic_1, y = statistic_2, fill = probability)
  ) +
    ggplot2::geom_raster() +
    ggplot2::geom_hline(yintercept = 0) +
    ggplot2::geom_vline(xintercept = 0) +
    ggplot2::coord_fixed() +
    ggplot2::labs(
      x = expression(delta[x]), y = expression(delta[y]),
      title = "Map of rejection probabilities", subtitle = subtitle,
      fill = "rejection\nprobability"
    )

  if (!is.null(filename)) {
    ggplot2::ggsave(filename = filename, plot = figure)
  } else {
    suppressWarnings(print(figure))
  }
  if (return_fig) figure else invisible(NULL)
}
