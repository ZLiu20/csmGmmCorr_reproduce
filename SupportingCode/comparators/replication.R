run_repfdr_replication_R1 <- function(z, t_value, q_value) {
  if (!requireNamespace("repfdr", quietly = TRUE)) stop("The repfdr package is required")
  z <- as.matrix(z)
  K <- ncol(z)
  binned <- repfdr::ztobins(z, n.association.status = 3, n.bins = 120,
                            plot.diagnostics = FALSE, trim.z = TRUE,
                            trim.z.upper = 8, trim.z.lower = -8)
  configurations <- as.matrix(repfdr::hconfigs(K, n.association.status = 3))
  nonzero <- rowSums(configurations != 0)
  concordant <- rowSums(configurations != 0) == abs(rowSums(sign(configurations)))
  alternative_rows <- which(nonzero > t_value & concordant)
  fitted <- repfdr::repfdr(zscores.binned = binned, non.null = "user.defined",
                           non.null.rows = alternative_rows, clusters.verb = FALSE)
  fdr <- clamp_probability_R1(fitted$mat[, "fdr"])
  list(rejected = fdr <= q_value, score = fdr, rankingScoreType = "local_fdr",
       tuning = paste0("user.defined replication,t=", t_value))
}

run_dirpc_replication_R1 <- function(z, t_value, q_value) {
  z <- as.matrix(z)
  K <- ncol(z)
  rank <- as.integer(t_value) + 1L
  pplus <- stats::pnorm(-z)
  pminus <- stats::pnorm(z)
  kth <- function(P) apply(P, 1L, function(v) sort(v, partial = rank)[rank])
  plus <- pmin(1, (K - t_value) * kth(pplus))
  minus <- pmin(1, (K - t_value) * kth(pminus))
  p <- pmin(1, 2 * pmin(plus, minus))
  adjusted <- stats::p.adjust(p, "BH")
  list(rejected = adjusted <= q_value, score = p, rankingScoreType = "p_value",
       tuning = paste0("directional PC,t=", t_value))
}
