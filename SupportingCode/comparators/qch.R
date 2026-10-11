run_qch_pleiotropy_R1 <- function(z, t_value, q_value, threads = 1L) {
  if (!requireNamespace("qch", quietly = TRUE)) stop("The qch package is required")
  p <- two_sided_p_R1(z)
  K <- ncol(p)
  at_least <- as.integer(t_value) + 1L
  null_config <- qch::GetHconfig(K)
  alt_config <- qch::GetH1AtLeast(null_config, AtLeast = at_least)
  test_name <- paste0("AtLeast_", at_least)
  fitted <- qch::qch.fit(pValMat = as.data.frame(p), Hconfig = null_config,
                         copula = "gaussian", threads_nb = as.integer(threads), plotting = FALSE)
  tested <- qch::qch.test(res.qch.fit = fitted, Hconfig = null_config,
                          Hconfig.H1 = alt_config, Alpha = q_value,
                          threads_nb = as.integer(threads))
  score <- clamp_probability_R1(tested$lFDR[[test_name]])
  list(rejected = as.logical(tested$Rejection[[test_name]]), score = score,
       rankingScoreType = "local_fdr",
       tuning = paste0("AtLeast=", at_least, ";copula=gaussian"))
}
