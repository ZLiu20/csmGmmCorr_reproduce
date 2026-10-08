# Data generation for the csmGmmCorr simulation study.
# Follows Sun, McCaw and Lin (csmGmm) and the pilot design of this study: n = 1000 subjects,
# 10^5 SNPs in sets of 1000, MAF 0.3, logistic outcomes with intercept -1, outcomes
# correlated through bindata::rmvbin-type latent normals with AR(1) correlation 0.2,
# Wald statistics from logistic regression truncated at +-8.1, and methods fitted and
# evaluated on the same data. Revised to match the paper:
#   * configurations are drawn iid from the two-group model: support size s with
#     probability lambda_s (s = 0..K, all positive), support uniform among the
#     choose(K, s) sets, signs uniform (conditional symmetry);
#   * mediation uses one common binary outcome with all mediators entering jointly.

# ---- design points ----
# Each scenario has its own grid (as in Sun, McCaw and Lin, whose pleiotropy, replication
# and mediation figures use different effect sizes): an effect-size axis (effect_01..05)
# and a signal-proportion axis (tau1_*, at the largest effect size). The grids were fixed
# before the structure-B runs by one rule: the power of the methods spans roughly 0.2-0.9
# (from pilot runs), so that no design is uninformative (power near 0 or 1).
#   pl  (pleiotropy):  beta = 0.38, 0.43, 0.48, 0.53, 0.58;  tau1 = 0.01-0.04 at 0.58
#   rp  (replication): beta = 0.48, 0.54, 0.60, 0.66, 0.72;  tau1 = 0.01-0.04 at 0.72
#       (replication needs equal signs at every non-null coordinate and has lower
#        power than pleiotropy at the same effect size)
#   med (mediation):   alpha in [0.00, 0.18] + 0.01 m, beta in [0.14, 0.24] + 0.01 m,
#       m = 8, 10, 12, 14, 16 (the form of Sun et al.); tau1 = 0.02-0.08 at m = 12
# study "main": latent AR(1) correlation rho = 0.2 (as in Sun et al.), intercept b0 = -1;
# study "rho": the effect_05 setting of each scenario with rho = 0, 0.2, ..., 0.8
#   for the outcomes (pl, rp) and the mediator errors (med);
# study "sym": the main pl / rp designs with intercept b0 = 0 (outcome prevalence 1/2),
#   for which the statistics of positive and negative effects have the same distribution
#   up to sign (with b0 = -1 positive effects give about 10-15 percent larger |z|, which
#   the conditionally symmetric model does not allow).
.grid <- list(
  pl = list(effect = c(0.38, 0.43, 0.48, 0.53, 0.58), tau1 = c(0.01, 0.02, 0.03, 0.04)),
  rp = list(effect = c(0.48, 0.54, 0.60, 0.66, 0.72), tau1 = c(0.01, 0.02, 0.03, 0.04)),
  med = list(m = c(8, 10, 12, 14, 16), m_tau1 = 12, tau1 = c(0.02, 0.04, 0.06, 0.08)))
# Extension of the effect-size axis to stronger signals (designs 10-12, appended so that the
# indices and seeds of designs 1-9 are unchanged): effect_06..08.
.grid_ext <- list(pl = c(0.63, 0.68, 0.73), rp = c(0.78, 0.84, 0.90), med = c(18, 20, 22))
design_table <- function(study = "main", scen = "pl") {
  if (!scen %in% names(.grid)) stop("unknown scenario: ", scen)
  g <- .grid[[scen]]
  if (scen == "med") {
    m <- c(g$m, rep(g$m_tau1, length(g$tau1)))
    main <- data.frame(effect = NA_real_, x = 0.01 * m, a_min = 0.01 * m, a_max = 0.18 + 0.01 * m,
                       b_min = 0.14 + 0.01 * m, b_max = 0.24 + 0.01 * m)
  } else {
    main <- data.frame(effect = c(g$effect, rep(max(g$effect), length(g$tau1))),
                       x = c(g$effect, rep(max(g$effect), length(g$tau1))),
                       a_min = NA, a_max = NA, b_min = NA, b_max = NA)
  }
  main <- cbind(data.frame(design = c(sprintf("effect_%02d", 1:5), sprintf("tau1_%03d", round(100 * g$tau1))),
                           axis = c(rep("effect", 5), rep("tau1", length(g$tau1))), stringsAsFactors = FALSE),
                main, data.frame(tau1 = c(rep(0, 5), g$tau1), rho = 0.2, b0 = -1, scen = scen, stringsAsFactors = FALSE))
  ge <- .grid_ext[[scen]]
  ext <- if (scen == "med") {
    data.frame(effect = NA_real_, x = 0.01 * ge, a_min = 0.01 * ge, a_max = 0.18 + 0.01 * ge,
               b_min = 0.14 + 0.01 * ge, b_max = 0.24 + 0.01 * ge)
  } else data.frame(effect = ge, x = ge, a_min = NA, a_max = NA, b_min = NA, b_max = NA)
  ext <- cbind(data.frame(design = sprintf("effect_%02d", 6:8), axis = "effect", stringsAsFactors = FALSE), ext,
               data.frame(tau1 = 0, rho = 0.2, b0 = -1, scen = scen, stringsAsFactors = FALSE))
  main$x[main$axis == "tau1"] <- main$tau1[main$axis == "tau1"]   # x: value on the design's axis
  main <- rbind(main, ext[names(main)])
  main$outcome <- "logit"; main$ld <- 0; main$corr <- "ar1"
  if (study %in% c("main", "split")) return(main)
  if (study == "sym") { main$b0 <- 0; return(main) }   # intercept 0: symmetric statistics (all scenarios)
  # Robustness studies (same designs and seeds layout, different data-generating process):
  #   quant: quantitative outcome(s), linear model with unit-variance errors and the effects
  #          scaled by sqrt(0.197) = 0.444, which gives the noncentrality of the logistic
  #          Wald statistic at prevalence logit^-1(-1); statistics are symmetric in sign;
  #   ld   : genotypes within a set of 1,000 SNPs in linkage disequilibrium (latent AR(1)
  #          correlation 0.5 between adjacent SNPs on each haplotype), so that statistics
  #          are dependent across features;
  #   exch : exchangeable latent correlation 0.4 between outcomes (mediator errors)
  #          instead of AR(1) with rho = 0.2.
  if (study == "quant") { main$outcome <- "linear"; return(main) }
  if (study %in% c("ld", "ldloc", "ldu")) { main$ld <- 0.5; return(main) }   # ldloc: same data, locus-level FDP
  if (study == "exch") { main$corr <- "exch"; main$rho <- 0.4; return(main) }
  if (study == "rho") {
    rho <- c(0, 0.2, 0.4, 0.6, 0.8)
    r <- main[rep(5L, length(rho)), ]; r$x <- rho
    r$design <- sprintf("rho_%03d", round(100 * rho)); r$axis <- "rho"; r$rho <- rho; r$tau1 <- 0
    rownames(r) <- NULL
    return(r)
  }
  stop("unknown study: ", study)
}

# Support-size probabilities lambda_0..lambda_K (all positive). Sizes 0-4 as in the pilot design;
# sizes 5..K continue the size-4 value with a halving tail.
size_probs <- function(K, tau1 = 0) {
  if (tau1 > 0) {
    p11 <- 3 * tau1^2
    head <- c(3 * tau1, 3 * p11, p11 / 2, p11 / 4)                 # sizes 1..4 (signal-proportion designs)
  } else {
    head <- c(0.04, 0.0012, 0.0006, 0.0003)                         # sizes 1..4 (effect-size designs)
  }
  lam <- head[seq_len(min(K, 4))]
  if (K > 4) lam <- c(lam, head[4] * 0.5^(seq_len(K - 4)))
  c(1 - sum(lam), lam)
}

# Fraction of size-s supports that contain at least one complete pair (2g-1, 2g).
complete_pair_fraction <- function(K, s) {
  G <- K / 2
  if (s < 2) return(0)
  none <- if (s <= G) choose(G, s) * 2^s else 0
  1 - none / choose(K, s)
}

# Mediation: rescale the sizes that can be alternatives (s > t) so that the
# expected proportion of mediation alternatives equals the pilot target
# (0.001 at effect points, 0.025 * tau1 at signal-proportion points).
size_probs_mediation <- function(K, t, tau1 = 0) {
  lam <- size_probs(K, tau1)
  target <- if (tau1 > 0) 0.025 * tau1 else 0.001
  s <- 0:K
  frac <- vapply(s, function(x) if (x > t) complete_pair_fraction(K, x) else 0, 0)
  base <- sum(lam * frac)
  elig <- s > t
  lam[elig] <- lam[elig] * target / base
  lam[1] <- 1 - sum(lam[-1])
  if (lam[1] <= 0) stop("mediation calibration gives a non-positive null probability")
  lam
}

# iid configurations: rows of {-1, 0, 1}^K.
draw_configs <- function(N, K, lam) {
  s <- sample.int(K + 1L, N, replace = TRUE, prob = lam) - 1L
  H <- matrix(0L, N, K)
  for (i in which(s > 0L)) {
    a <- sample.int(K, s[i])
    H[i, a] <- sample(c(-1L, 1L), s[i], replace = TRUE)
  }
  H
}

ar1 <- function(K, rho) rho^abs(outer(seq_len(K), seq_len(K), "-"))
corr_mat <- function(K, rho, corr = "ar1") {
  if (corr == "exch") { R <- matrix(rho, K, K); diag(R) <- 1; R } else ar1(K, rho)
}
LIN_SCALE <- sqrt(plogis(-1) * (1 - plogis(-1)))   # 0.444: linear effect with the logistic noncentrality

# Genotypes (n x m, 0/1/2) with MAF maf. ld = 0: independent SNPs. ld > 0: each haplotype is a
# thresholded latent Gaussian AR(1) sequence along the m SNPs with lag-one correlation ld.
gen_geno <- function(n, m, maf, ld = 0) {
  if (ld <= 0) return(matrix(rbinom(n * m, 2L, maf), n, m))
  thr <- qnorm(maf); s <- sqrt(1 - ld^2)
  hap <- function() {
    U <- matrix(rnorm(n * m), n, m)
    for (j in 2:m) U[, j] <- ld * U[, j - 1L] + s * U[, j]
    (U < thr) * 1L
  }
  hap() + hap()
}

# Score statistic sqrt(n) * cor(g, y) of y ~ 1 + g for every column g of G (linear model).
lin_score_z <- function(G, y) {
  Gc <- sweep(G, 2L, colMeans(G)); yc <- y - mean(y)
  sqrt(length(y)) * colSums(Gc * yc) / sqrt(colSums(Gc^2) * sum(yc^2))
}

# Batched least squares with t statistics: one model per column j, Y[, j] ~ X_1[, j] + ... (Xs as in
# logit_wald_batch); returns the p x m matrix of coefficient / standard error.
lin_wald_batch <- function(Xs, Y) {
  n <- nrow(Y); m <- ncol(Y); p <- length(Xs)
  z <- matrix(0, p, m)
  for (j in seq_len(m)) {
    X <- vapply(Xs, function(x) x[, j], numeric(n))
    f <- stats::lm.fit(X, Y[, j])
    if (f$rank < p || any(f$qr$pivot != seq_len(p))) next
    s2 <- sum(f$residuals^2) / (n - p)
    R <- qr.R(f$qr); Ri <- backsolve(R, diag(p))
    z[, j] <- f$coefficients / sqrt(s2 * rowSums(Ri^2))
  }
  z
}

Z_MAX <- 8.1    # statistics truncated at +-8.1, as in Sun et al.

# Batched logistic regression with Wald statistics: one model per column j,
# Y[, j] ~ X_1[, j] + ... + X_p[, j] (Xs: list of n x m covariate matrices, the first
# being the intercept). Newton-Raphson; returns the p x m coefficients, standard
# errors and Wald statistics coef / se (as summary(glm(...))$coefficients[, 3]).
logit_wald_batch <- function(Xs, Y, maxit = 25L, tol = 1e-8) {
  n <- nrow(Y); m <- ncol(Y); p <- length(Xs)
  th <- matrix(0, p, m); th[1, ] <- qlogis(pmin(pmax(colMeans(Y), 1e-6), 1 - 1e-6))
  lin <- function(th) { e <- matrix(0, n, m); for (k in seq_len(p)) e <- e + Xs[[k]] * rep(th[k, ], each = n); e }
  info <- function(w) {
    I <- array(0, c(p, p, m))
    for (u in seq_len(p)) for (v in u:p) { val <- colSums(w * Xs[[u]] * Xs[[v]]); I[u, v, ] <- val; I[v, u, ] <- val }
    I
  }
  inv <- function(I) {                       # p x p x m inverses
    out <- array(NA_real_, dim(I))
    if (p == 2L) {
      det <- I[1, 1, ] * I[2, 2, ] - I[1, 2, ]^2
      out[1, 1, ] <- I[2, 2, ] / det; out[2, 2, ] <- I[1, 1, ] / det; out[1, 2, ] <- out[2, 1, ] <- -I[1, 2, ] / det
    } else for (j in seq_len(m)) out[, , j] <- tryCatch(solve(I[, , j]), error = function(e) matrix(NA_real_, p, p))
    out
  }
  for (it in seq_len(maxit)) {
    pr <- plogis(lin(th)); w <- pr * (1 - pr); r <- Y - pr
    S <- t(matrix(vapply(Xs, function(x) colSums(x * r), numeric(m)), ncol = p))   # p x m
    Iv <- inv(info(w))
    D <- matrix(0, p, m)
    for (u in seq_len(p)) for (v in seq_len(p)) D[u, ] <- D[u, ] + Iv[u, v, ] * S[v, ]
    D[!is.finite(D)] <- 0
    th <- th + D
    if (max(abs(D)) < tol) break
  }
  pr <- plogis(lin(th)); Iv <- inv(info(pr * (1 - pr)))
  se <- sqrt(pmax(t(vapply(seq_len(p), function(u) Iv[u, u, ], numeric(m))), 0))
  if (m == 1L) se <- matrix(se, p, 1L)
  z <- th / se
  z[!is.finite(z)] <- 0
  list(coef = th, se = se, z = z)
}

# Wald z of y ~ 1 + g for every column g of G (genotypes 0/1/2) and a common outcome y.
# With a three-level covariate the likelihood depends on the data only through the
# counts n_c and events e_c at g = c (c = 0, 1, 2), so Newton-Raphson runs on three
# cells per SNP; the estimate, standard error and Wald z equal those of
# glm(y ~ g, family = binomial) on the individual data.
logit_wald_geno <- function(G, y, maxit = 50L, tol = 1e-10) {
  m <- ncol(G)
  nc <- rbind(colSums(G == 0), colSums(G == 1), colSums(G == 2))       # 3 x m
  ec <- rbind(colSums((G == 0) * y), colSums((G == 1) * y), colSums((G == 2) * y))
  g <- c(0, 1, 2)
  a <- rep(qlogis(min(max(mean(y), 1e-6), 1 - 1e-6)), m); b <- rep(0, m)
  for (it in seq_len(maxit)) {
    p <- plogis(rbind(a, a + b, a + 2 * b)); w <- nc * p * (1 - p); r <- ec - nc * p
    i11 <- colSums(w); i12 <- colSums(w * g); i22 <- colSums(w * g^2)
    s1 <- colSums(r); s2 <- colSums(r * g)
    det <- i11 * i22 - i12^2
    da <- (i22 * s1 - i12 * s2) / det; db <- (i11 * s2 - i12 * s1) / det
    da[!is.finite(da)] <- 0; db[!is.finite(db)] <- 0
    a <- a + da; b <- b + db
    if (max(abs(c(da, db))) < tol) break
  }
  p <- plogis(rbind(a, a + b, a + 2 * b)); w <- nc * p * (1 - p)
  i11 <- colSums(w); i12 <- colSums(w * g); i22 <- colSums(w * g^2)
  z <- b / sqrt(i11 / (i11 * i22 - i12^2))
  z[!is.finite(z)] <- 0
  z
}

# ---- pleiotropy / replication data (shared; fit once, scored for t = 1, 2) ----
# As in Sun et al. (Fig. 8 code): within a set of set_size SNPs, outcome k of
# subject i has eta_ik = b0 + sum_j G_ij beta_jk (all SNPs of the set jointly); the K
# outcomes of a subject are generated as in bindata::rmvbin(margprob, sigma = AR1(rho)),
# i.e. thresholded latent normals with AR(1) correlation rho; the statistic of SNP j for
# outcome k is the Wald z of y_k ~ 1 + g_j.
gen_pleio <- function(N, K, lam, effect, n = 1000L, maf = 0.3, b0 = -1, rho = 0.2,
                      set_size = 1000L, outcome = "logit", ld = 0, corr = "ar1") {
  H <- draw_configs(N, K, lam)
  B <- H * effect
  Lc <- chol(corr_mat(K, rho, corr))
  Z <- matrix(0, N, K)
  for (st in seq(1L, N, by = set_size)) {
    idx <- st:min(N, st + set_size - 1L); m <- length(idx)
    G <- gen_geno(n, m, maf, ld)
    Lat <- matrix(rnorm(n * K), n, K) %*% Lc
    if (outcome == "linear") {
      Y <- G %*% (LIN_SCALE * B[idx, , drop = FALSE]) + Lat                # n x K, correlated errors
      for (k in seq_len(K)) Z[idx, k] <- lin_score_z(G, Y[, k])
    } else {
      eta <- G %*% B[idx, , drop = FALSE] + b0                       # n x K
      Y <- (Lat <= qnorm(plogis(eta))) * 1
      for (k in seq_len(K)) Z[idx, k] <- logit_wald_geno(G, Y[, k])
    }
  }
  list(Z = pmax(pmin(Z, Z_MAX), -Z_MAX), H = H, B = B)
}

# ---- mediation data (common outcome, joint mediators) ----
# For SNP j: M_jg = alpha_jg G_j + e_jg (e AR(1) across mediators g),
#            logit P(Y_j = 1) = b0 + sum_g beta_jg M_jg (no direct effect);
# statistics: alpha_g -> score statistic of G in the linear model M_g ~ G (as in Sun et al.);
#             beta_g -> Wald z of M_g in the joint logistic model Y ~ G + M_1 + ... + M_G
#             (Sun et al. use the Wald z of M in Y ~ G + M with one outcome per mediator).
# Columns are ordered (alpha_1, beta_1, alpha_2, beta_2, ...).
gen_mediation <- function(N, K, lam, a_rng, b_rng, n = 1000L, maf = 0.3, b0 = -1,
                          rho = 0.2, set_size = 1000L, outcome = "logit", ld = 0, corr = "ar1") {
  Gm <- K / 2
  H <- draw_configs(N, K, lam)
  odd <- seq(1L, K, by = 2L); even <- odd + 1L
  B <- matrix(0, N, K)
  B[, odd] <- H[, odd] * matrix(runif(N * Gm, a_rng[1], a_rng[2]), N, Gm)
  B[, even] <- H[, even] * matrix(runif(N * Gm, b_rng[1], b_rng[2]), N, Gm)
  Lc <- chol(corr_mat(Gm, rho, corr))
  Z <- matrix(0, N, K)
  for (st in seq(1L, N, by = set_size)) {
    idx <- st:min(N, st + set_size - 1L); m <- length(idx)
    G <- gen_geno(n, m, maf, ld)
    E <- matrix(rnorm(n * m * Gm), n * m, Gm) %*% Lc
    M <- lapply(seq_len(Gm), function(g)
      G * rep(B[idx, odd[g]], each = n) + matrix(E[, g], n, m))
    eta <- matrix(if (outcome == "linear") 0 else b0, n, m)
    bsc <- if (outcome == "linear") LIN_SCALE else 1
    for (g in seq_len(Gm)) eta <- eta + M[[g]] * rep(bsc * B[idx, even[g]], each = n)
    Gc <- sweep(G, 2L, colMeans(G))
    for (g in seq_len(Gm)) {
      Mc <- sweep(M[[g]], 2L, colMeans(M[[g]]))
      Z[idx, odd[g]] <- sqrt(n) * colSums(Gc * Mc) / sqrt(colSums(Gc^2) * colSums(Mc^2))
    }
    if (outcome == "linear") {
      Y <- eta + matrix(rnorm(n * m), n, m)
      wz <- lin_wald_batch(c(list(matrix(1, n, m), G), M), Y)
    } else {
      Y <- matrix(rbinom(n * m, 1L, plogis(eta)), n, m)
      wz <- logit_wald_batch(c(list(matrix(1, n, m), G), M), Y)$z     # rows: 1, G, M_1..M_Gm
    }
    for (g in seq_len(Gm)) Z[idx, even[g]] <- wz[2L + g, ]
  }
  Z[!is.finite(Z)] <- 0
  list(Z = pmax(pmin(Z, Z_MAX), -Z_MAX), H = H, B = B)
}

# Seeds: one stream per (scenario, design, K, t, replicate, role).
sim_seed <- function(scen, design_id, K, t, rep, role) {
  # Replicates 1-200 keep the original seeds. Replicates above 200 use a wider
  # replicate field, which keeps every seed distinct from those of replicates 1-200.
  if (rep > 200L)
    return(as.integer(1000000000 + ((((scen * 10L + design_id) * 11L + K) * 3L + t) * 1000L + rep) * 2L + role) %% .Machine$integer.max)
  as.integer(((((scen * 10L + design_id) * 11L + K) * 3L + t) * 100L + rep) * 2L + role) %% .Machine$integer.max
}

# One dataset per replicate; as in Sun et al., every method is fitted and evaluated
# on the same data (train and test below are the same object). Seed streams:
# pl 1 / 3 / 5, med 2 / 4, rp 7 / 9 / 11 (main / rho / sym).
make_dataset <- function(scen, d, K, t, rep, N = 100000L, study = "main") {
  # study "split": the training data are those of the main study (same design
  # table, same seeds), and the test data are an independent dataset of the same design
  # (role 2), as in the setting of main-text Theorem 3.
  split <- study == "split"; if (split) study <- "main"
  dt <- design_table(study, scen)[d, ]
  stream <- c(pl = 1L, med = 2L, rp = 7L)[[scen]] +
    switch(study, main = 0L, rho = 2L, sym = 4L, quant = 20L, ld = 40L, ldloc = 40L, ldu = 40L, exch = 60L)
  gen <- function(role) {
    set.seed(sim_seed(stream, d, K, t, rep, role))
    if (scen %in% c("pl", "rp")) {
      gen_pleio(N, K, size_probs(K, dt$tau1), dt$effect, rho = dt$rho, b0 = dt$b0,
                outcome = dt$outcome, ld = dt$ld, corr = dt$corr)
    } else {
      gen_mediation(N, K, size_probs_mediation(K, t, dt$tau1), c(dt$a_min, dt$a_max),
                    c(dt$b_min, dt$b_max), rho = dt$rho, b0 = dt$b0,
                    outcome = dt$outcome, ld = dt$ld, corr = dt$corr)
    }
  }
  dat <- gen(1L)
  list(train = dat, test = if (split) gen(2L) else dat, design = dt)
}

truth_alt <- function(H, target, t) {
  s <- rowSums(H != 0)
  switch(target,
    pleiotropy = s > t,
    replication = s > t & abs(rowSums(H)) == s,
    mediation = {
      pr <- matrix(seq_len(ncol(H)), ncol = 2L, byrow = TRUE)
      s > t & rowSums((H[, pr[, 1], drop = FALSE] != 0) & (H[, pr[, 2], drop = FALSE] != 0)) >= 1
    })
}

# Locus-level truth (LD study): feature j counts as a true discovery at the locus level
# if some feature of the same set of set_size consecutive SNPs within w positions of j is an
# alternative (lag-w latent correlation 0.5^w; w = 5 gives 0.03).
truth_alt_locus <- function(alt, w = 5L, set_size = 1000L) {
  N <- length(alt); out <- logical(N)
  for (st in seq(1L, N, by = set_size)) {
    idx <- st:min(N, st + set_size - 1L); a <- alt[idx]; m <- length(a)
    cs <- c(0L, cumsum(a)); lo <- pmax(1L, seq_len(m) - w); hi <- pmin(m, seq_len(m) + w)
    out[idx] <- (cs[hi + 1L] - cs[lo]) > 0L
  }
  out
}

# Locus-level truth by union of nearby effects: the effect configuration of feature j is
# taken as the signed union of the configurations of the features of its set within w positions
# (coordinate k nonzero if some such feature has a nonzero effect on k, with the sign of their sum),
# so that a null variant in linkage disequilibrium with, e.g., two single-trait variants of different
# traits counts as pleiotropic. The target is then evaluated on this union configuration.
union_config <- function(H, w = 2L, set_size = 1000L) {
  N <- nrow(H); U <- matrix(0L, N, ncol(H))
  for (st in seq(1L, N, by = set_size)) {
    idx <- st:min(N, st + set_size - 1L); m <- length(idx)
    for (k in seq_len(ncol(H))) {
      h <- H[idx, k]; cs <- c(0, cumsum(h)); ca <- c(0L, cumsum(h != 0))
      lo <- pmax(1L, seq_len(m) - w); hi <- pmin(m, seq_len(m) + w)
      any_nz <- (ca[hi + 1L] - ca[lo]) > 0L; sm <- cs[hi + 1L] - cs[lo]
      U[idx, k] <- ifelse(any_nz, ifelse(sm >= 0, 1L, -1L), 0L)
    }
  }
  U
}
