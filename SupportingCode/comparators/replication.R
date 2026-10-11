# Comparison methods for the replication target (more than t nonzero effects, all of one sign).
#
# repfdr (Heller and Yekutieli, 2014): empirical Bayes over the 3^K association-status vectors
#   (-1, 0, +1 per study) with study-wise binned densities and independent studies; the
#   non-null set is user defined as exactly the replication alternative (more than t nonzero
#   entries, all of the same sign). Features are rejected when the estimated Bayes FDR
#   (Fdr, the running mean of the sorted local fdr) is at most q. Package repfdr (CRAN).
# Directional partial conjunction: one-sided p-values p+ = 1 - Phi(z), p- = Phi(z);
#   P+ = min(1, (K - t) p+_(t+1)), P- likewise, P = min(1, 2 min(P+, P-)), then BH.
#   It tests "more than t effects in one direction" and ignores effects of the opposite
#   sign, so it is NOT valid for the replication null, which also contains discordant
#   configurations; it is included to show the consequence.
run_repfdr_replication_R1 <- function(z, t_value, q_value) {
  if (!requireNamespace("repfdr", quietly = TRUE)) stop("The repfdr package is required.")
  z <- as.matrix(z); K <- ncol(z)
  zb <- repfdr::ztobins(z, n.association.status = 3, n.bins = 120, plot.diagnostics = FALSE,
                        trim.z = TRUE, trim.z.upper = 8, trim.z.lower = -8)
  H <- as.matrix(repfdr::hconfigs(K, n.association.status = 3))
  nz <- rowSums(H != 0); conc <- abs(rowSums(sign(H))) == nz
  rows <- which(nz > t_value & conc)
  fit <- repfdr::repfdr(zb$pdf.binned.z, zb$binned.z.mat, non.null = "user.defined",
                        non.null.rows = rows, clusters.verbose = FALSE)
  fdr <- clamp_probability_R1(fit$mat[, "fdr"]); Fdr <- clamp_probability_R1(fit$mat[, "Fdr"])
  list(rejected = as.logical(Fdr <= q_value), score = fdr, rankingScoreType = "local_fdr",
       tuning = paste0("user.defined replication, t=", t_value))
}
run_dirpc_replication_R1 <- function(z, t_value, q_value) {
  z <- as.matrix(z); K <- ncol(z); k <- as.integer(t_value) + 1L
  pp <- stats::pnorm(-z); pm <- stats::pnorm(z)
  kth <- function(P) apply(P, 1L, function(v) sort(v, partial = k)[k])
  Pp <- pmin(1, (K - t_value) * kth(pp)); Pm <- pmin(1, (K - t_value) * kth(pm))
  P <- pmin(1, 2 * pmin(Pp, Pm))
  adj <- stats::p.adjust(P, "BH")
  list(rejected = as.logical(adj <= q_value), score = P, rankingScoreType = "p_value",
       tuning = paste0("directional PC, t=", t_value))
}
