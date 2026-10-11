.simulation_grid <- list(
  pl = list(effect = c(0.38, 0.43, 0.48, 0.53, 0.58), tau = c(0.01, 0.02, 0.03, 0.04)),
  rp = list(effect = c(0.48, 0.54, 0.60, 0.66, 0.72), tau = c(0.01, 0.02, 0.03, 0.04)),
  med = list(m = c(8, 10, 12, 14, 16), m_tau = 12,
             tau = c(0.02, 0.04, 0.06, 0.08))
)

.strong_effects <- list(pl = c(0.63, 0.68, 0.73),
                        rp = c(0.78, 0.84, 0.90),
                        med = c(18, 20, 22))

design_table <- function(study = "main", scen = "pl") {
  if (!scen %in% names(.simulation_grid)) stop("Unknown scenario: ", scen)
  grid <- .simulation_grid[[scen]]

  if (scen == "med") {
    m <- c(grid$m, rep(grid$m_tau, length(grid$tau)))
    core <- data.frame(
      effect = NA_real_,
      x = 0.01 * m,
      a_min = 0.01 * m,
      a_max = 0.18 + 0.01 * m,
      b_min = 0.14 + 0.01 * m,
      b_max = 0.24 + 0.01 * m
    )
  } else {
    e <- c(grid$effect, rep(max(grid$effect), length(grid$tau)))
    core <- data.frame(effect = e, x = e, a_min = NA_real_, a_max = NA_real_,
                       b_min = NA_real_, b_max = NA_real_)
  }

  labels <- c(sprintf("effect_%02d", 1:5), sprintf("tau1_%03d", round(100 * grid$tau)))
  ans <- cbind(
    data.frame(design = labels,
               axis = c(rep("effect", 5L), rep("tau1", length(grid$tau))),
               stringsAsFactors = FALSE),
    core,
    data.frame(tau1 = c(rep(0, 5L), grid$tau), rho = 0.2, b0 = -1,
               scen = scen, stringsAsFactors = FALSE)
  )

  stronger <- .strong_effects[[scen]]
  if (scen == "med") {
    extra <- data.frame(
      effect = NA_real_, x = 0.01 * stronger,
      a_min = 0.01 * stronger, a_max = 0.18 + 0.01 * stronger,
      b_min = 0.14 + 0.01 * stronger, b_max = 0.24 + 0.01 * stronger
    )
  } else {
    extra <- data.frame(effect = stronger, x = stronger, a_min = NA_real_, a_max = NA_real_,
                        b_min = NA_real_, b_max = NA_real_)
  }
  extra <- cbind(
    data.frame(design = sprintf("effect_%02d", 6:8), axis = "effect",
               stringsAsFactors = FALSE),
    extra,
    data.frame(tau1 = 0, rho = 0.2, b0 = -1, scen = scen,
               stringsAsFactors = FALSE)
  )
  ans <- rbind(ans, extra)
  ans$x[ans$axis == "tau1"] <- ans$tau1[ans$axis == "tau1"]

  ans$outcome <- "logit"
  ans$ld <- 0
  ans$corr <- "ar1"
  if (study %in% c("main", "split")) return(ans)
  if (study == "sym") {
    ans$b0 <- 0
    return(ans)
  }
  if (study == "quant") {
    ans$outcome <- "linear"
    return(ans)
  }
  if (study %in% c("ld", "ldloc", "ldu")) {
    ans$ld <- 0.5
    return(ans)
  }
  if (study == "exch") {
    ans$corr <- "exch"
    ans$rho <- 0.4
    return(ans)
  }
  if (study == "rho") {
    rho <- seq(0, 0.8, by = 0.2)
    ans <- ans[rep(5L, length(rho)), , drop = FALSE]
    ans$design <- sprintf("rho_%03d", round(100 * rho))
    ans$axis <- "rho"
    ans$x <- rho
    ans$rho <- rho
    ans$tau1 <- 0
    rownames(ans) <- NULL
    return(ans)
  }
  stop("Unknown simulation study: ", study)
}

size_probs <- function(K, tau1 = 0) {
  if (tau1 > 0) {
    p11 <- 3 * tau1^2
    leading <- c(3 * tau1, 3 * p11, p11 / 2, p11 / 4)
  } else {
    leading <- c(0.04, 0.0012, 0.0006, 0.0003)
  }
  signal <- leading[seq_len(min(K, 4L))]
  if (K > 4L) signal <- c(signal, leading[4L] * 0.5^seq_len(K - 4L))
  c(1 - sum(signal), signal)
}

complete_pair_fraction <- function(K, support_size) {
  groups <- K / 2
  if (support_size < 2) return(0)
  without_pair <- if (support_size <= groups) choose(groups, support_size) * 2^support_size else 0
  1 - without_pair / choose(K, support_size)
}

size_probs_mediation <- function(K, t, tau1 = 0) {
  probs <- size_probs(K, tau1)
  target <- if (tau1 > 0) 0.025 + tau1 else 0.001
  sizes <- 0:K
  eligible <- sizes > t
  pair_share <- vapply(sizes, function(s) {
    if (s > t) complete_pair_fraction(K, s) else 0
  }, numeric(1))
  present <- sum(probs * pair_share)
  probs[eligible] <- probs[eligible] * target / present
  probs[1L] <- 1 - sum(probs[-1L])
  if (probs[1L] <= 0) stop("Mediation calibration produced an invalid null mass")
  probs
}

draw_configs <- function(N, K, probabilities) {
  support_size <- sample.int(K + 1L, N, replace = TRUE, prob = probabilities) - 1L
  H <- matrix(0L, N, K)
  for (i in which(support_size > 0L)) {
    active <- sample.int(K, support_size[i], replace = FALSE)
    H[i, active] <- sample(c(-1L, 1L), support_size[i], replace = TRUE)
  }
  H
}

.correlation_matrix <- function(K, rho, type = "ar1") {
  if (type == "exch") {
    R <- matrix(rho, K, K)
    diag(R) <- 1
    return(R)
  }
  rho^abs(outer(seq_len(K), seq_len(K), "-"))
}

.linear_scale <- sqrt(plogis(-1) * (1 - plogis(-1)))
.z_limit <- 8.1

gen_genotypes <- function(n, m, maf = 0.3, ld = 0) {
  if (ld <= 0) return(matrix(rbinom(n * m, 2L, maf), n, m))
  latent <- matrix(0, n, m)
  threshold <- qnorm(maf)
  innovation <- sqrt(1 - ld^2)
  latent[, 1L] <- rnorm(n)
  if (m > 1L) for (j in 2:m) latent[, j] <- ld * latent[, j - 1L] + innovation * rnorm(n)
  haplotype <- function() (latent < threshold) * 1L
  haplotype() + haplotype()
}

linear_score_z <- function(G, y) {
  centered_G <- sweep(G, 2L, colMeans(G))
  centered_y <- y - mean(y)
  numerator <- colSums(centered_G * centered_y)
  sqrt(length(y)) * numerator / sqrt(colSums(centered_G^2) * sum(centered_y^2))
}

linear_wald_batch <- function(Xs, Y) {
  n <- nrow(Y)
  m <- ncol(Y)
  p <- length(Xs)
  z <- matrix(0, p, m)
  for (j in seq_len(m)) {
    X <- vapply(Xs, function(x) x[, j], numeric(n))
    fit <- qr(X)
    if (fit$rank < p) next
    coef <- qr.coef(fit, Y[, j])
    residual <- qr.resid(fit, Y[, j])
    R <- qr.R(fit)
    inverse <- chol2inv(R)
    sigma2 <- sum(residual^2) / (n - p)
    se <- sqrt(diag(inverse) * sigma2)
    z[, j] <- coef / se
  }
  z[!is.finite(z)] <- 0
  z
}

logit_wald_batch <- function(Xs, Y, maxit = 25L, tol = 1e-8) {
  n <- nrow(Y)
  m <- ncol(Y)
  p <- length(Xs)
  theta <- matrix(0, p, m)
  theta[1L, ] <- qlogis(pmin(pmax(colMeans(Y), 1e-6), 1 - 1e-6))

  information <- function(weights) {
    result <- array(0, c(p, p, m))
    for (u in seq_len(p)) for (v in u:p) {
      value <- colSums(weights * Xs[[u]] * Xs[[v]])
      result[u, v, ] <- value
      result[v, u, ] <- value
    }
    result
  }
  invert_pages <- function(A) {
    result <- array(NA_real_, dim(A))
    if (p == 2L) {
      determinant <- A[1, 1, ] * A[2, 2, ] - A[1, 2, ]^2
      result[1, 1, ] <- A[2, 2, ] / determinant
      result[2, 2, ] <- A[1, 1, ] / determinant
      result[1, 2, ] <- result[2, 1, ] <- -A[1, 2, ] / determinant
    } else {
      for (j in seq_len(m)) result[, , j] <- tryCatch(solve(A[, , j]), error = function(e) matrix(NA_real_, p, p))
    }
    result
  }

  for (iteration in seq_len(maxit)) {
    eta <- matrix(0, n, m)
    for (u in seq_len(p)) eta <- eta + sweep(Xs[[u]], 2L, theta[u, ], `*`)
    fitted <- plogis(eta)
    weights <- fitted * (1 - fitted)
    score <- vapply(Xs, function(x) colSums(x * (Y - fitted)), numeric(m))
    score <- t(score)
    inverse <- invert_pages(information(weights))
    step <- matrix(0, p, m)
    for (j in seq_len(m)) step[, j] <- inverse[, , j] %*% score[, j]
    step[!is.finite(step)] <- 0
    theta <- theta + step
    if (max(abs(step)) < tol) break
  }

  eta <- matrix(0, n, m)
  for (u in seq_len(p)) eta <- eta + sweep(Xs[[u]], 2L, theta[u, ], `*`)
  inverse <- invert_pages(information(plogis(eta) * (1 - plogis(eta))))
  se <- matrix(NA_real_, p, m)
  for (j in seq_len(m)) se[, j] <- sqrt(diag(inverse[, , j]))
  z <- theta / se
  z[!is.finite(z)] <- 0
  list(coef = theta, se = se, z = z)
}

logit_wald_genotype <- function(G, y, maxit = 50L, tol = 1e-10) {
  m <- ncol(G)
  counts <- rbind(colSums(G == 0L), colSums(G == 1L), colSums(G == 2L))
  events <- rbind(colSums((G == 0L) * y), colSums((G == 1L) * y), colSums((G == 2L) * y))
  intercept <- rep(qlogis(pmin(pmax(mean(y), 1e-6), 1 - 1e-6)), m)
  slope <- numeric(m)
  genotype <- c(0, 1, 2)
  for (iteration in seq_len(maxit)) {
    probability <- plogis(outer(genotype, slope, `*`) + rep(intercept, each = 3L))
    weights <- counts * probability * (1 - probability)
    residual <- events - counts * probability
    i11 <- colSums(weights)
    i12 <- colSums(weights * genotype)
    i22 <- colSums(weights * genotype^2)
    s1 <- colSums(residual)
    s2 <- colSums(residual * genotype)
    determinant <- i11 * i22 - i12^2
    da <- (i22 * s1 - i12 * s2) / determinant
    db <- (i11 * s2 - i12 * s1) / determinant
    da[!is.finite(da)] <- 0
    db[!is.finite(db)] <- 0
    intercept <- intercept + da
    slope <- slope + db
    if (max(abs(c(da, db))) < tol) break
  }
  probability <- plogis(outer(genotype, slope, `*`) + rep(intercept, each = 3L))
  weights <- counts * probability * (1 - probability)
  i11 <- colSums(weights)
  i12 <- colSums(weights * genotype)
  i22 <- colSums(weights * genotype^2)
  z <- slope / sqrt(i11 / (i11 * i22 - i12^2))
  z[!is.finite(z)] <- 0
  z
}

gen_pleiotropy <- function(N, K, probabilities, effect, n = 1000L, maf = 0.3,
                           b0 = -1, rho = 0.2, set_size = 1000L,
                           outcome = "logit", ld = 0, corr = "ar1") {
  H <- draw_configs(N, K, probabilities)
  beta <- H * effect
  root <- chol(.correlation_matrix(K, rho, corr))
  Z <- matrix(0, N, K)
  for (start in seq.int(1L, N, by = set_size)) {
    idx <- start:min(N, start + set_size - 1L)
    m <- length(idx)
    G <- gen_genotypes(n, m, maf, ld)
    latent <- matrix(rnorm(n * K), n, K) %*% root
    if (outcome == "linear") {
      Y <- G %*% (.linear_scale * beta[idx, , drop = FALSE]) + latent
      for (k in seq_len(K)) Z[idx, k] <- linear_score_z(G, Y[, k])
    } else {
      eta <- G %*% beta[idx, , drop = FALSE] + b0
      Y <- (latent <= qnorm(plogis(eta))) * 1L
      for (k in seq_len(K)) Z[idx, k] <- logit_wald_genotype(G, Y[, k])
    }
  }
  list(Z = pmax(pmin(Z, .z_limit), -.z_limit), H = H, B = beta)
}

gen_mediation <- function(N, K, probabilities, alpha_range, beta_range,
                          n = 1000L, maf = 0.3, b0 = -1, rho = 0.2,
                          set_size = 1000L, outcome = "logit", ld = 0,
                          corr = "ar1") {
  mediators <- K / 2L
  H <- draw_configs(N, K, probabilities)
  alpha_idx <- seq.int(1L, K, by = 2L)
  beta_idx <- alpha_idx + 1L
  B <- matrix(0, N, K)
  B[, alpha_idx] <- H[, alpha_idx] * matrix(runif(N * mediators, alpha_range[1L], alpha_range[2L]), N, mediators)
  B[, beta_idx] <- H[, beta_idx] * matrix(runif(N * mediators, beta_range[1L], beta_range[2L]), N, mediators)
  root <- chol(.correlation_matrix(mediators, rho, corr))
  Z <- matrix(0, N, K)

  for (start in seq.int(1L, N, by = set_size)) {
    idx <- start:min(N, start + set_size - 1L)
    m <- length(idx)
    G <- gen_genotypes(n, m, maf, ld)
    errors <- matrix(rnorm(n * m * mediators), n * m, mediators) %*% root
    mediator_values <- vector("list", mediators)
    eta <- matrix(if (outcome == "linear") 0 else b0, n, m)
    beta_scale <- if (outcome == "linear") .linear_scale else 1

    for (g in seq_len(mediators)) {
      mediator_values[[g]] <- G * rep(B[idx, alpha_idx[g]], each = n) + matrix(errors[, g], n, m)
      eta <- eta + mediator_values[[g]] * rep(beta_scale * B[idx, beta_idx[g]], each = n)
      centered_G <- sweep(G, 2L, colMeans(G))
      centered_M <- sweep(mediator_values[[g]], 2L, colMeans(mediator_values[[g]]))
      Z[idx, alpha_idx[g]] <- colSums(centered_G * centered_M) /
        sqrt(colSums(centered_G^2) * colSums(centered_M^2)) * sqrt(n - 2)
    }

    if (outcome == "linear") {
      Y <- eta + matrix(rnorm(n * m), n, m)
      fit <- linear_wald_batch(c(list(matrix(1, n, m), G), mediator_values), Y)
    } else {
      Y <- matrix(rbinom(n * m, 1L, plogis(eta)), n, m)
      fit <- logit_wald_batch(c(list(matrix(1, n, m), G), mediator_values), Y)$z
    }
    for (g in seq_len(mediators)) Z[idx, beta_idx[g]] <- fit[2L + g, ]
  }
  Z[!is.finite(Z)] <- 0
  list(Z = pmax(pmin(Z, .z_limit), -.z_limit), H = H, B = B)
}

sim_seed <- function(stream, design, K, t, replicate, role) {
  key <- ((((stream * 10L + design) * 11L + K) * 3L + t) * 1000L + replicate) * 2L + role
  if (replicate > 200L) {
    return(as.integer((1000000000 + key) %% .Machine$integer.max))
  }
  as.integer(key %% .Machine$integer.max)
}

make_dataset <- function(scen, d, K, t, replicate, N = 100000L, study = "main") {
  split_sample <- identical(study, "split")
  effective_study <- if (split_sample) "main" else study
  design <- design_table(effective_study, scen)[d, , drop = FALSE]
  stream <- c(pl = 1L, med = 2L, rp = 7L)[[scen]] +
    c(main = 0L, rho = 2L, sym = 4L, quant = 20L, ld = 40L,
      ldloc = 40L, ldu = 40L, exch = 60L)[[effective_study]]

  generate <- function(role) {
    set.seed(sim_seed(stream, d, K, t, replicate, role))
    if (scen %in% c("pl", "rp")) {
      gen_pleiotropy(N, K, size_probs(K, design$tau1), design$effect,
                     rho = design$rho, b0 = design$b0, outcome = design$outcome,
                     ld = design$ld, corr = design$corr)
    } else {
      gen_mediation(N, K, size_probs_mediation(K, t, design$tau1),
                    c(design$a_min, design$a_max), c(design$b_min, design$b_max),
                    rho = design$rho, b0 = design$b0, outcome = design$outcome,
                    ld = design$ld, corr = design$corr)
    }
  }

  train <- generate(1L)
  list(train = train, test = if (split_sample) generate(2L) else train, design = design)
}

truth_alt <- function(H, target, t) {
  active <- rowSums(H != 0)
  switch(target,
         pleiotropy = active > t,
         replication = active > t & abs(rowSums(H)) == active,
         mediation = {
           pair_index <- matrix(seq_len(ncol(H)), ncol = 2L, byrow = TRUE)
           complete <- rowSums((H[, pair_index[, 1L], drop = FALSE] != 0) &
                                 (H[, pair_index[, 2L], drop = FALSE] != 0))
           active > t & complete >= 1L
         },
         stop("Unknown target: ", target))
}

truth_alt_locus <- function(alt, window = 5L, set_size = 1000L) {
  out <- logical(length(alt))
  for (start in seq.int(1L, length(alt), by = set_size)) {
    idx <- start:min(length(alt), start + set_size - 1L)
    cumulative <- c(0L, cumsum(alt[idx]))
    lo <- pmax(1L, seq_along(idx) - window)
    hi <- pmin(length(idx), seq_along(idx) + window)
    out[idx] <- cumulative[hi + 1L] - cumulative[lo] > 0L
  }
  out
}

union_config <- function(H, window = 2L, set_size = 1000L) {
  out <- matrix(0L, nrow(H), ncol(H))
  for (start in seq.int(1L, nrow(H), by = set_size)) {
    idx <- start:min(nrow(H), start + set_size - 1L)
    for (k in seq_len(ncol(H))) {
      values <- H[idx, k]
      signed_sum <- c(0L, cumsum(values))
      nonzero_sum <- c(0L, cumsum(values != 0L))
      lo <- pmax(1L, seq_along(idx) - window)
      hi <- pmin(length(idx), seq_along(idx) + window)
      has_signal <- nonzero_sum[hi + 1L] - nonzero_sum[lo] > 0L
      direction <- signed_sum[hi + 1L] - signed_sum[lo]
      out[idx, k] <- ifelse(has_signal, ifelse(direction >= 0, 1L, -1L), 0L)
    }
  }
  out
}
