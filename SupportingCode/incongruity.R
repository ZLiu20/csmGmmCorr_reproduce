.order_representation <- function(z, target, t) {
  if (target == "replication") {
    rank <- t + 1L
    ordered <- function(M) apply(M, 1L, function(v) sort(v, decreasing = TRUE)[rank])
    positive <- ordered(z) >= ordered(-z)
    list(group = ifelse(positive, "+", "-"),
         evidence = z * ifelse(positive, 1, -1))
  } else {
    signs <- as.vector(z > 0) * 2^(rep(seq_len(ncol(z)) - 1L, each = nrow(z)))
    group <- rowSums(matrix(signs, nrow(z), ncol(z)))
    zero <- rowSums(z == 0) > 0
    group <- as.character(group)
    group[zero] <- paste0("z", which(zero))
    list(group = group, evidence = abs(z))
  }
}

incon_decision <- function(z, rejected, target, t) {
  rejected <- as.logical(rejected)
  if (!any(rejected) || all(rejected)) return(0L)
  representation <- .order_representation(z, target, t)
  W <- representation$evidence
  K <- ncol(W)
  flagged <- logical(nrow(W))
  for (group in unique(representation$group[rejected])) {
    yes <- which(rejected & representation$group == group)
    no <- which(!rejected & representation$group == group)
    if (!length(no)) next
    lower <- apply(W[yes, , drop = FALSE], 2L, min)
    no <- no[rowSums(W[no, , drop = FALSE] >= rep(lower, each = length(no))) == K]
    if (!length(no)) next
    candidates <- W[no, , drop = FALSE]
    for (j in yes) {
      boundary <- matrix(W[j, ], nrow(candidates), K, byrow = TRUE)
      dominates <- rowSums(candidates >= boundary) == K & rowSums(candidates > boundary) > 0L
      if (any(dominates)) flagged[no[dominates]] <- TRUE
    }
  }
  sum(flagged)
}

incon_score <- function(z, score, target, t,
                        max_n = as.integer(Sys.getenv("INCON_MAX", "5000"))) {
  if (is.null(score) || all(is.na(score))) return(NA_integer_)
  keep <- which(is.finite(score) & score < 0.99)
  if (length(keep) <= 1L) return(0L)
  keep <- keep[order(score[keep])][seq_len(min(length(keep), max_n))]
  representation <- .order_representation(z[keep, , drop = FALSE], target, t)
  W <- representation$evidence
  s <- score[keep]
  violations <- 0L
  for (group in unique(representation$group)) {
    members <- which(representation$group == group)
    if (length(members) < 2L) next
    members <- members[do.call(order, c(list(s[members]),
                                       lapply(seq_len(ncol(W)), function(k) -W[members, k])))]
    ordered_W <- W[members, , drop = FALSE]
    for (i in 2:length(members)) {
      previous <- ordered_W[seq_len(i - 1L), , drop = FALSE]
      current <- matrix(ordered_W[i, ], i - 1L, ncol(W), byrow = TRUE)
      if (any(rowSums(previous <= current) == ncol(W) & rowSums(previous < current) > 0L)) {
        violations <- violations + 1L
      }
    }
  }
  violations
}
