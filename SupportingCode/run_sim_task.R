# One task of a simulation array: replicates rep_from..rep_to of one (design, K, t) cell.
#   run_sim_task("main_ours_pl_pen", aID)   (array names: SupportingCode/sim_tasks.R)
# Every replicate regenerates its data from fixed seeds (sim_data.R, sim_seed), fits the
# methods of the array and writes one CSV with one row per (target, t, method):
#   <folder>/output/<study>/<type>/<design>_K<K>_t<t>_rep<rep>.csv
# Replicates whose file exists are skipped, so a task can be resubmitted.
# Environment variables: SIM_N (features per data set, default 100000; other values write
# to output_smoke/ and are for testing only), MAX_ITER (ECM iterations, default 1500),
# LSB_DJOB_NUMPROC (threads; set by LSF from -n), SIM_MAX_REPS (testing: replicates per task).
run_sim_task <- function(array, aID, root = here::here()) {
  A <- sim_array(array)
  if (aID < 1L || aID > nrow(A$tasks)) stop(sprintf("task %d outside 1-%d", aID, nrow(A$tasks)))
  tk <- A$tasks[aID, ]
  study <- A$study; type0 <- A$type; scen <- A$scen
  d <- tk$design; K <- tk$K; t_data <- tk$t
  cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "1"))
  Q <- 0.1
  N <- as.integer(Sys.getenv("SIM_N", "100000"))
  MAXIT <- as.integer(Sys.getenv("MAX_ITER", "1500"))
  dt <- design_table(study, scen)[d, ]
  out_dir <- file.path(root, A$folder, if (N != 100000L) "output_smoke" else "output", study, type0)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  if (grepl("^ours|_csm$", type0)) suppressMessages(library(csmGmmCorr))
  message(sprintf("array %s, task %d: design %d (%s), K = %d, t = %d, replicates %d-%d, N = %d, cores = %d",
                  array, aID, d, dt$design, K, t_data, tk$rep_from, tk$rep_to, N, cores))

  run_rep <- function(r) {
    out <- file.path(out_dir, sprintf("%s_K%02d_t%d_rep%03d.csv", dt$design, K, t_data, r))
    if (file.exists(out)) { message("exists: ", out); return(invisible()) }
    t0 <- proc.time()[["elapsed"]]
    dat <- make_dataset(scen, d, K, t_data, r, N, study = study)
    gen_min <- (proc.time()[["elapsed"]] - t0) / 60
    te <- dat$test; tr <- dat$train
    targets <- switch(scen, pl = "pleiotropy", rp = "replication", med = "mediation")
    ts <- if (scen == "med") t_data else 1:2
    base <- data.frame(study = study, scen = scen, design = dt$design, axis = dt$axis, x = dt$x,
                       effect = dt$effect, tau1 = dt$tau1, a_min = dt$a_min, rho = dt$rho, K = K, rep = r,
                       stringsAsFactors = FALSE)
    row <- function(target, t, method, rej, extra = list(), score = NULL) {
      alt <- truth_alt(te$H, target, t)
      o <- cbind(base, data.frame(target = target, t = t, method = method, nAlt = sum(alt), R = sum(rej),
                 FP = sum(rej & !alt), FDP = if (sum(rej)) mean(!alt[rej]) else 0,
                 power = if (sum(alt)) sum(rej & alt) / sum(alt) else NA,
                 incon_dec = incon_decision(te$Z, rej, target, t),
                 incon_score = incon_score(te$Z, score, target, t), stringsAsFactors = FALSE))
      if (isTRUE(dt$ld > 0)) {   # LD study: locus-level FDP and union-of-effects truth
        al <- truth_alt_locus(alt, 5L)
        o$FP_loc <- sum(rej & !al); o$FDP_loc <- if (sum(rej)) mean(!al[rej]) else 0
        for (w in c(1L, 2L, 5L)) {
          au <- truth_alt(union_config(te$H, w), target, t)
          o[[paste0("FDP_union", w)]] <- if (sum(rej)) mean(!au[rej]) else 0
          o[[paste0("power_union", w)]] <- if (sum(au)) sum(rej & au) / sum(au) else NA
        }
      }
      if (length(extra)) o <- cbind(o, as.data.frame(extra, stringsAsFactors = FALSE))
      o
    }
    res <- NULL
    fit_or_fail <- function(Z, sp, ctl)
      tryCatch(csm_fit(Z, sp, starts = c("null", "shared"), control = ctl),
               error = function(e) list(converged = FALSE, message = conditionMessage(e)))

    if (grepl("^ours", type0)) {
      # Region: magnitudes in [0.5, 15], kappa = 0.05. Magnitude components: structure B (2 for
      # supports of size one, 1 otherwise); A = 1/1 and C = 2/2 for the sensitivity arrays.
      struct <- if (grepl("_[ABC]$", type0)) sub("^.*_", "", type0) else "B"
      nc <- switch(struct, A = c(1L, 1L), B = c(2L, 1L), C = c(2L, 2L))
      pen <- grepl("_pen", type0)
      lam <- 30; alp <- 10; pl_tag <- ""
      if (grepl("_L[0-9]+A[0-9]+$", type0)) {
        lam <- as.numeric(sub("^.*_L([0-9]+)A([0-9]+)$", "\\1", type0))
        alp <- as.numeric(sub("^.*_L([0-9]+)A([0-9]+)$", "\\2", type0))
        pl_tag <- sprintf("_L%gA%g", lam, alp)
      }
      ctl <- if (pen) csm_control(cores = cores, max_iter = MAXIT, penalty = lam, weight_penalty = alp)
             else csm_control(cores = cores, max_iter = MAXIT, penalty = 0, weight_penalty = 0)
      v <- paste0(if (pen) "penalized" else "ML", if (struct != "B") paste0("_", struct) else "", pl_tag)
      if (grepl("_oracle$", type0)) {
        # Sigma fixed at the projected correlation of the statistics of the null features.
        R0 <- csmGmmCorr:::.project_cor(stats::cor(tr$Z[rowSums(tr$H != 0) == 0, , drop = FALSE]), 0.05)
        ctl <- csm_control(cores = cores, max_iter = MAXIT, penalty = lam, weight_penalty = alp, fix_sigma = R0)
        v <- "penalized_nullSigma"
      }
      sp <- csm_spec(K, M = csm_structure(K, nc[1], nc[2]), b_mu = 15)
      tm <- system.time(f <- fit_or_fail(tr$Z, sp, ctl))[["elapsed"]]
      message(sprintf("replicate %d: %s fit %.1f min, converged %s", r, v, tm / 60, f$converged))
      info <- list(converged = f$converged,
                   iterations = if (f$converged) sum(f$runs$iterations) else NA,
                   fit_min = tm / 60, gen_min = gen_min,
                   Sigma_minEig = if (f$converged) min(eigen(f$params$Sigma, TRUE, TRUE)$values) else NA,
                   structure = struct,
                   loglik = if (f$converged) as.numeric(logLik(f)) else NA,
                   df = if (f$converged) attr(logLik(f), "df") else NA)
      for (tg in targets) for (t in ts) {
        if (!f$converged) { res <- rbind(res, row(tg, t, v, rep(FALSE, N), c(info, list(flagged = NA)))); next }
        lf <- csm_lfdr(f, te$Z, tg, t, cores = cores)
        su <- csm_stepup(lf, Q)
        ex <- csm_extrapolation(f, te$Z, tg, t, su$rejected, su$cutoff)
        res <- rbind(res, row(tg, t, v, su$rejected, c(info, list(flagged = sum(ex$flag, na.rm = TRUE))), score = lf))
        ad <- csm_adjust(te$Z, lf, tg, t, Q)
        res <- rbind(res, row(tg, t, paste0(v, "_adjusted"), ad$rejected, c(info, list(flagged = NA))))
      }
    } else if (grepl("_csm$", type0)) {
      # csmGmm (Sun et al.): the same model with Sigma fixed, structure B.
      sp <- csm_spec(K, M = csm_structure(K, 2L, 1L), b_mu = 15)
      R_plug <- csmGmmCorr:::.project_cor(stats::cor(tr$Z), sp$kappa)
      for (v in c("csmGmm", "csmGmm_cor", "csmGmm_pen")) {
        ctl <- switch(v,
          csmGmm     = csm_control(cores = cores, max_iter = MAXIT, penalty = 0, weight_penalty = 0, fix_sigma = TRUE),
          csmGmm_cor = csm_control(cores = cores, max_iter = MAXIT, penalty = 0, weight_penalty = 0, fix_sigma = R_plug),
          csmGmm_pen = csm_control(cores = cores, max_iter = MAXIT, fix_sigma = TRUE))
        tm <- system.time(f <- fit_or_fail(tr$Z, sp, ctl))[["elapsed"]]
        message(sprintf("replicate %d: %s fit %.1f min, converged %s", r, v, tm / 60, f$converged))
        info <- list(converged = f$converged, iterations = if (f$converged) sum(f$runs$iterations) else NA,
                     fit_min = tm / 60, gen_min = gen_min, structure = "B")
        for (tg in targets) for (t in ts) {
          if (!f$converged) { res <- rbind(res, row(tg, t, v, rep(FALSE, N), info)); next }
          lf <- csm_lfdr(f, te$Z, tg, t, cores = cores)
          su <- csm_stepup(lf, Q)
          res <- rbind(res, row(tg, t, v, su$rejected, info, score = lf))
        }
      }
    } else {
      comp <- switch(type0,
        comp_pl  = list(target = "pleiotropy", fns = list(QCH = function(Z, t) run_qch_pleiotropy_R1(Z, t, Q, threads = cores),
                                                          adaFilter = function(Z, t) run_adafilter_pleiotropy_R1(Z, t, Q))),
        comp_rp  = list(target = "replication", fns = list(repfdr = function(Z, t) run_repfdr_replication_R1(Z, t, Q),
                                                           DirectionalPC = function(Z, t) run_dirpc_replication_R1(Z, t, Q))),
        comp_med = list(target = "mediation", fns = list(DEI = function(Z, t) run_dei_mediation_R1(Z, t, Q),
                                                         MediationTest = function(Z, t) run_mediation_test_R1(Z, t, Q),
                                                         PCPathSafeguard = function(Z, t) run_mediation_safeguard_R1(Z, t, Q))),
        stop("unknown method type: ", type0))
      for (t in ts) for (m in names(comp$fns)) {
        o <- tryCatch(comp$fns[[m]](te$Z, t), error = function(e) e)
        res <- rbind(res, if (inherits(o, "error")) cbind(row(comp$target, t, m, rep(FALSE, N)), error = conditionMessage(o))
                          else cbind(row(comp$target, t, m, o$rejected, score = o$score), error = ""))
      }
    }
    res$total_min <- (proc.time()[["elapsed"]] - t0) / 60
    utils::write.csv(res, out, row.names = FALSE)
    message("done: ", out)
  }
  reps <- tk$rep_from:tk$rep_to
  nmax <- as.integer(Sys.getenv("SIM_MAX_REPS", "0"))   # testing only: first nmax replicates of the task
  if (nmax > 0L) reps <- utils::head(reps, nmax)
  for (r in reps) run_rep(r)
  invisible(NULL)
}
