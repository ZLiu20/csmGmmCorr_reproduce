# QCH (Mary-Huard et al.): pleiotropy, H1 = at least t + 1 non-null coordinates.
run_qch_pleiotropy_R1 <- function(z, t_value, q_value, threads = 1L) {
  if (!requireNamespace("qch", quietly = TRUE)) {
    stop("The qch package is required for the qch_copula comparator.")
  }
  p_matrix <- two_sided_p_R1(z)
  K <- ncol(p_matrix)
  at_least <- as.integer(t_value) + 1L
  h_config <- qch::GetHconfig(K)
  h1_config <- qch::GetH1AtLeast(h_config, AtLeast = at_least)
  test_name <- paste0("AtLeast_", at_least)
  fit <- qch::qch.fit(
    pValMat = as.data.frame(p_matrix),
    Hconfig = h_config,
    copula = "gaussian",
    threads_nb = as.integer(threads),
    plotting = FALSE
  )
  tested <- qch::qch.test(
    res.qch.fit = fit,
    Hconfig = h_config,
    Hconfig.H1 = h1_config,
    Alpha = q_value,
    threads_nb = as.integer(threads)
  )
  score <- clamp_probability_R1(tested$lFDR[[test_name]])
  rejected <- as.logical(tested$Rejection[[test_name]])
  list(
    rejected = rejected,
    score = score,
    rankingScoreType = "local_fdr",
    tuning = paste0("AtLeast=", at_least, ";copula=gaussian")
  )
}

