FLOOR_DESIGN <- list(
  K = 4L, probs = c(0.93, 0.05, 0.012, 0.005, 0.003),
  mu = c(2.5, 3, 2.5, 3), rho = 0.3,
  N = 20000L, J = 20000L, t = 1L, q = 0.1
)

floor_configs <- function(K) {
  H <- as.matrix(expand.grid(rep(list(c(-1L, 0L, 1L)), K)))
  dimnames(H) <- NULL
  list(H = H, S = rowSums(H != 0L))
}

floor_simulate <- function(n, probs, mu, Sigma) {
  K <- length(mu)
  support <- sample.int(K + 1L, n, replace = TRUE, prob = probs) - 1L
  H <- matrix(0L, n, K)
  for (i in which(support > 0L)) {
    active <- sample.int(K, support[i])
    H[i, active] <- sample(c(-1L, 1L), support[i], replace = TRUE)
  }
  noise <- matrix(stats::rnorm(n * K), n, K) %*% chol(Sigma)
  list(Z = H * rep(mu, each = n) + noise, H = H)
}

floor_logdens <- function(Z, configs, omega, mu, Sigma) {
  U <- chol(Sigma)
  K <- ncol(Z)
  weight <- omega[configs$S + 1L] /
    choose(K, configs$S) / pmax(1, 2^configs$S)
  vapply(seq_len(nrow(configs$H)), function(i) {
    residual <- backsolve(U, t(Z) - configs$H[i, ] * mu, transpose = TRUE)
    log(weight[i]) - 0.5 * colSums(residual^2) -
      sum(log(diag(U))) - 0.5 * K * log(2 * pi)
  }, numeric(nrow(Z)))
}

floor_weights <- function(n_s, lower) {
  if (lower <= 0) return(n_s / sum(n_s))
  objective <- function(lambda) sum(pmax(lower, n_s / lambda)) - 1
  root <- stats::uniroot(objective, c(1e-12, 10 * sum(n_s) + 1), tol = 1e-14)$root
  out <- pmax(lower, n_s / root)
  out / sum(out)
}

floor_fit <- function(Z, a_omega, max_iter = 2000L, tol = 1e-9, a_mu = 0.5) {
  K <- ncol(Z)
  n <- nrow(Z)
  cf <- floor_configs(K)
  omega <- floor_weights(c(0.8, 0.1, 0.05, 0.03, 0.02)[seq_len(K + 1L)] * n, a_omega)
  mu <- rep(2, K)
  Sigma <- stats::cor(Z)
  old <- -Inf
  for (iteration in seq_len(max_iter)) {
    logd <- floor_logdens(Z, cf, omega, mu, Sigma)
    row_max <- apply(logd, 1L, max)
    prob <- exp(logd - row_max)
    total <- rowSums(prob)
    ll <- mean(row_max + log(total))
    resp <- prob / total
    count <- colSums(resp)
    weighted_z <- crossprod(resp, Z)
    omega <- floor_weights(vapply(0:K, function(s) sum(count[cf$S == s]), 0), a_omega)
    precision <- solve(Sigma)
    A <- matrix(0, K, K)
    b <- numeric(K)
    for (r in seq_len(nrow(cf$H))) {
      h <- cf$H[r, ]
      A <- A + count[r] * (outer(h, h) * precision)
      b <- b + h * as.vector(precision %*% weighted_z[r, ])
    }
    mu <- pmax(a_mu, as.vector(solve(A + diag(1e-10, K), b)))
    mean_states <- cf$H * rep(mu, each = nrow(cf$H))
    scatter <- crossprod(Z)
    for (r in seq_len(nrow(cf$H))) {
      m <- mean_states[r, ]
      v <- weighted_z[r, ]
      scatter <- scatter - tcrossprod(m, v) - tcrossprod(v, m) + count[r] * tcrossprod(m)
    }
    Sigma <- stats::cov2cor((scatter + t(scatter)) / (2 * n))
    if (abs(ll - old) < tol) break
    old <- ll
  }
  list(omega = omega, mu = mu, Sigma = Sigma, loglik = ll, iterations = iteration)
}

floor_lfdr <- function(Z, omega, mu, Sigma, target, t = 1L) {
  cf <- floor_configs(ncol(Z))
  logd <- floor_logdens(Z, cf, omega, mu, Sigma)
  p <- exp(logd - apply(logd, 1L, max))
  null <- !csm_truth(cf$H, target, t)
  rowSums(p[, null, drop = FALSE]) / rowSums(p)
}

floor_stepup <- function(lfdr, q) csm_stepup(lfdr, q)$rejected

floor_study_rep <- function(rep, floors = c(1 / 64, 1e-4)) {
  d <- FLOOR_DESIGN
  Sigma <- matrix(d$rho, d$K, d$K)
  diag(Sigma) <- 1
  set.seed(20261007 + rep)
  training <- floor_simulate(d$N, d$probs, d$mu, Sigma)
  testing <- floor_simulate(d$J, d$probs, d$mu, Sigma)
  fits <- c(list(oracle = list(omega = d$probs, mu = d$mu, Sigma = Sigma)),
            stats::setNames(lapply(floors, function(a) floor_fit(training$Z, a)), sprintf("floor_%g", floors)))
  out <- NULL
  for (method in names(fits)) for (target in c("pleiotropy", "replication", "mediation")) {
    f <- fits[[method]]
    lfdr <- floor_lfdr(testing$Z, f$omega, f$mu, f$Sigma, target, d$t)
    reject <- floor_stepup(lfdr, d$q)
    alt <- csm_truth(testing$H, target, d$t)
    out <- rbind(out, data.frame(
      replicate = rep, fit = method, target = target, R = sum(reject),
      FDP = if (any(reject)) mean(!alt[reject]) else 0,
      power = if (any(alt)) sum(reject & alt) / sum(alt) else NA_real_,
      omega2to4 = sum(f$omega[3:5]),
      iterations = if (is.null(f$iterations)) NA_integer_ else f$iterations
    ))
  }
  out
}
