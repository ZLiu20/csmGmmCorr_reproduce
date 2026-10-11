# Incongruity counts for every method (base R only).
#
# Order of evidence (as in the paper and csmGmmCorr::csm_adjust): feature j is at least as
# strong as feature i when both are in the same class and w_j >= w_i coordinatewise.
#   pleiotropy, mediation: class = sign orthant of z (a zero coordinate: own class), w = |z|
#     (the order used by Sun, McCaw and Lin for incongruous results);
#   replication: class = directional branch (+ if the (t+1)-th largest z is at least the
#     (t+1)-th largest -z, else -), w = z on + and -z on -.
# Two counts:
#   incon_dec   : decision level. Number of features that are NOT rejected although they
#                 dominate (>= in every coordinate, > in at least one) a rejected feature of
#                 the same class. Applies to any method; zero for the monotone-adjusted rule.
#   incon_score : score level, Sun et al.'s count. Among features with score < 0.99 (the
#                 local fdr for our method and QCH; the ranking p-value for the other
#                 methods), restricted to the INCON_MAX (default 5000) smallest scores,
#                 the number of features whose score is larger than that of a feature they
#                 dominate in the same class.

.incon_repr <- function(Z, target, t) {
  if (target == "replication") {
    k <- t + 1L
    srt <- function(M) apply(M, 1L, function(v) sort(v, decreasing = TRUE)[k])
    plus <- srt(Z) >= srt(-Z)
    list(group = ifelse(plus, "+", "-"), w = Z * ifelse(plus, 1, -1))
  } else {
    g <- as.vector((Z > 0) %*% 2^(seq_len(ncol(Z)) - 1L))
    zero <- rowSums(Z == 0) > 0
    g <- as.character(g); g[zero] <- paste0("z", which(zero))
    list(group = g, w = abs(Z))
  }
}

incon_decision <- function(Z, rejected, target, t) {
  rejected <- as.logical(rejected)
  if (!any(rejected) || all(rejected)) return(0L)
  rp <- .incon_repr(Z, target, t); W <- rp$w; K <- ncol(W)
  flagged <- logical(nrow(W))
  for (gr in unique(rp$group[rejected])) {
    rj <- which(rejected & rp$group == gr); nr <- which(!rejected & rp$group == gr)
    if (!length(nr)) next
    lo <- apply(W[rj, , drop = FALSE], 2L, min)                   # can only dominate if >= this
    nr <- nr[rowSums(W[nr, , drop = FALSE] >= rep(lo, each = length(nr))) == K]
    if (!length(nr)) next
    Wn <- W[nr, , drop = FALSE]
    for (j in rj) {
      wj <- rep(W[j, ], each = length(nr))
      dom <- rowSums(Wn >= wj) == K & rowSums(Wn > wj) > 0
      if (any(dom)) flagged[nr[dom]] <- TRUE
    }
  }
  sum(flagged)
}

incon_score <- function(Z, score, target, t, max_n = as.integer(Sys.getenv("INCON_MAX", "5000"))) {
  if (is.null(score) || all(is.na(score))) return(NA_integer_)
  keep <- which(is.finite(score) & score < 0.99)
  if (length(keep) <= 1L) return(0L)
  keep <- keep[order(score[keep])][seq_len(min(length(keep), max_n))]
  rp <- .incon_repr(Z[keep, , drop = FALSE], target, t); W <- rp$w; K <- ncol(W); s <- score[keep]
  n_bad <- 0L
  for (gr in unique(rp$group)) {
    m <- which(rp$group == gr)
    if (length(m) < 2L) next
    m <- m[do.call(order, c(list(s[m]), lapply(seq_len(K), function(k) -W[m, k])))]   # Sun's ordering
    Wm <- W[m, , drop = FALSE]
    for (i in 2:length(m)) {
      prev <- Wm[seq_len(i - 1L), , drop = FALSE]; wi <- rep(Wm[i, ], each = i - 1L)
      if (any(rowSums(prev <= wi) == K & rowSums(prev < wi) > 0)) n_bad <- n_bad + 1L
    }
  }
  n_bad
}
