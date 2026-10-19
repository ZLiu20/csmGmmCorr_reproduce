# Supplementary Section S4.4.2 (Table S6): effect of a lower bound a_omega on the fitted
# support-size probabilities. Sparse design: K = 4, t = 1, q = 0.10, support-size probabilities
# (0.93, 0.05, 0.012, 0.005, 0.003), magnitudes (2.5, 3, 2.5, 3), equicorrelation 0.3,
# N = J = 20,000 independent training and test features per replicate.
# Fitted model: the exchangeable restriction (Supplementary Section S1.1), i.e. magnitudes
# mu_k shared by all supports and the probability omega_s of support size s spread uniformly
# over the supports and signs of that size, fitted by EM with omega_s >= a_omega for every s.
# Written independently of csmGmmCorr (which has no weight floor).
FLOOR_DESIGN <- list(K = 4L, probs = c(0.93, 0.05, 0.012, 0.005, 0.003), mu = c(2.5, 3, 2.5, 3), rho = 0.3,
                     N = 20000L, J = 20000L, t = 1L, q = 0.1)

floor_configs <- function(K) {
  H <- as.matrix(expand.grid(rep(list(c(-1L, 0L, 1L)), K))); dimnames(H) <- NULL
  list(H = H, S = rowSums(H != 0))
}

floor_simulate <- function(n, probs, mu, Sigma) {
  K <- length(mu)
  s <- sample.int(K + 1L, n, replace = TRUE, prob = probs) - 1L
  H <- matrix(0L, n, K)
  for (i in which(s > 0L)) { a <- sample.int(K, s[i]); H[i, a] <- sample(c(-1L, 1L), s[i], replace = TRUE) }
  Z <- H * rep(mu, each = n) + matrix(stats::rnorm(n * K), n, K) %*% chol(Sigma)
  list(Z = Z, H = H)
}

# Log of pi_h * phi(z; h * mu, Sigma) for all features (rows) and configurations (columns).
floor_logdens <- function(Z, cf, omega, mu, Sigma) {
  K <- ncol(Z); U <- chol(Sigma)
  w <- omega[cf$S + 1L] / choose(K, cf$S) / 2^cf$S
  sapply(seq_len(nrow(cf$H)), function(r) {
    D <- t(Z) - cf$H[r, ] * mu
    R <- backsolve(U, D, transpose = TRUE)
    log(w[r]) - 0.5 * colSums(R^2) - sum(log(diag(U))) - 0.5 * K * log(2 * pi)
  })
}

# Maximize sum_s n_s log(omega_s) subject to omega_s >= a and sum(omega) = 1:
# omega_s = max(a, n_s / lambda) with lambda chosen so that the omega sum to one.
floor_weights <- function(n_s, a) {
  if (a <= 0) return(n_s / sum(n_s))
  f <- function(l) sum(pmax(a, n_s / l)) - 1
  l <- stats::uniroot(f, c(1e-12, 10 * sum(n_s) + 1), tol = 1e-14)$root
  w <- pmax(a, n_s / l); w / sum(w)
}

floor_fit <- function(Z, a_omega, max_iter = 2000L, tol = 1e-9, a_mu = 0.5) {
  K <- ncol(Z); N <- nrow(Z); cf <- floor_configs(K)
  omega <- floor_weights(c(0.8, 0.1, 0.05, 0.03, 0.02)[seq_len(K + 1L)] * N, a_omega)
  mu <- rep(2, K); Sigma <- stats::cor(Z); ll_old <- -Inf
  for (it in seq_len(max_iter)) {
    L <- floor_logdens(Z, cf, omega, mu, Sigma)
    mx <- apply(L, 1, max); P <- exp(L - mx); tot <- rowSums(P)
    ll <- sum(mx + log(tot)) / N
    Rr <- P / tot                                       # responsibilities, N x configurations
    Rh <- colSums(Rr); Zh <- crossprod(Rr, Z)           # sums of r and of r z per configuration
    omega <- floor_weights(vapply(0:K, function(s) sum(Rh[cf$S == s]), 0), a_omega)
    Om <- solve(Sigma); A <- matrix(0, K, K); b <- numeric(K)
    for (r in seq_len(nrow(cf$H))) {
      h <- cf$H[r, ]
      A <- A + Rh[r] * (outer(h, h) * Om)
      b <- b + h * as.vector(Om %*% Zh[r, ])
    }
    mu <- pmax(a_mu, as.vector(solve(A, b)))
    M <- cf$H * rep(mu, each = nrow(cf$H))              # configuration means
    Sc <- (crossprod(Z) - crossprod(Zh, M) - crossprod(M, Zh) + crossprod(M * sqrt(Rh))) / N
    Sigma <- stats::cov2cor((Sc + t(Sc)) / 2)
    if (abs(ll - ll_old) < tol) break
    ll_old <- ll
  }
  list(omega = omega, mu = mu, Sigma = Sigma, loglik = ll, iterations = it)
}

floor_lfdr <- function(Z, omega, mu, Sigma, target, t = 1L) {
  cf <- floor_configs(ncol(Z)); L <- floor_logdens(Z, cf, omega, mu, Sigma)
  P <- exp(L - apply(L, 1, max))
  null <- !csm_truth(cf$H, target, t)
  rowSums(P[, null, drop = FALSE]) / rowSums(P)
}
csm_truth <- function(H, target, t) truth_alt(H, target, t)   # sim_data.R
floor_stepup <- function(lf, q) {
  o <- order(lf); k <- max(c(0L, which(cumsum(lf[o]) / seq_along(o) <= q)))
  rej <- logical(length(lf)); if (k > 0) rej[o[seq_len(k)]] <- TRUE; rej
}

floor_study_rep <- function(rep, floors = c(1 / 64, 1e-4)) {
  D <- FLOOR_DESIGN
  Sigma <- matrix(D$rho, D$K, D$K); diag(Sigma) <- 1
  set.seed(20261007 + rep)
  tr <- floor_simulate(D$N, D$probs, D$mu, Sigma); te <- floor_simulate(D$J, D$probs, D$mu, Sigma)
  fits <- c(list(oracle = list(omega = D$probs, mu = D$mu, Sigma = Sigma)),
            stats::setNames(lapply(floors, function(a) floor_fit(tr$Z, a)), sprintf("floor_%g", floors)))
  res <- NULL
  for (nm in names(fits)) for (tg in c("pleiotropy", "replication", "mediation")) {
    f <- fits[[nm]]
    lf <- floor_lfdr(te$Z, f$omega, f$mu, f$Sigma, tg, D$t)
    rej <- floor_stepup(lf, D$q); alt <- truth_alt(te$H, tg, D$t)
    res <- rbind(res, data.frame(rep = rep, fit = nm, target = tg, R = sum(rej),
                                 FDP = if (any(rej)) mean(!alt[rej]) else 0, power = sum(rej & alt) / sum(alt),
                                 omega2to4 = sum(f$omega[3:5]), iterations = if (is.null(f$iterations)) NA else f$iterations))
  }
  res
}
