run_sim_task <- function(array, aID, root = here::here()) {
  spec <- sim_array(array)
  aID <- as.integer(aID)
  if (aID < 1L || aID > nrow(spec$tasks)) {
    stop(sprintf("Task index %d is outside 1--%d", aID, nrow(spec$tasks)))
  }

  task <- spec$tasks[aID, , drop = FALSE]
  study <- spec$study
  method_type <- spec$type
  scenario <- spec$scen
  design_id <- task$design
  K <- task$K
  data_t <- task$t
  cores <- as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "1"))
  q <- 0.1
  N <- as.integer(Sys.getenv("SIM_N", "100000"))
  max_iter <- as.integer(Sys.getenv("MAX_ITER", "1500"))
  design <- design_table(study, scenario)[design_id, , drop = FALSE]

  result_dir <- file.path(root, if (N == 100000L) "output" else "output_smoke",
                          study, method_type)
  dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
  if (identical(method_type, "ours_rp_csm")) suppressMessages(library(csmGmmCorr))

  message(sprintf(
    "array=%s task=%d design=%d (%s) K=%d t=%d reps=%d-%d N=%d cores=%d",
    array, aID, design_id, design$design, K, data_t, task$rep_from, task$rep_to, N, cores
  ))

  evaluate_one <- function(replicate) {
    destination <- file.path(
      result_dir,
      sprintf("%s_K%02d_t%d_rep%03d.csv", design$design, K, data_t, replicate)
    )
    if (file.exists(destination)) {
      message("exists: ", destination)
      return(invisible(NULL))
    }

    started <- proc.time()[["elapsed"]]
    generated <- system.time(
      dataset <- make_dataset(scenario, design_id, K, data_t, replicate, N, study)
    )[["elapsed"]]
    train <- dataset$train
    test <- dataset$test
    targets <- c(pl = "pleiotropy", rp = "replication", med = "mediation")
    target <- targets[[scenario]]
    thresholds <- if (scenario == "med") data_t else 1:2

    common <- data.frame(
      study = study,
      scen = scenario,
      design = design$design,
      axis = design$axis,
      x = design$x,
      effect = design$effect,
      tau1 = design$tau1,
      a_min = design$a_min,
      rho = design$rho,
      K = K,
      rep = replicate,
      stringsAsFactors = FALSE
    )

    make_row <- function(target, t, method, rejected, extra = list(), score = NULL) {
      alt <- truth_alt(test$H, target, t)
      rejected <- as.logical(rejected)
      row <- cbind(
        common,
        data.frame(
          target = target,
          t = t,
          method = method,
          nAlt = sum(alt),
          R = sum(rejected),
          FP = sum(rejected & !alt),
          FDP = if (sum(rejected)) mean(!alt[rejected]) else 0,
          power = if (sum(alt)) sum(rejected & alt) / sum(alt) else NA_real_,
          incon_dec = incon_decision(test$Z, rejected, target, t),
          incon_score = incon_score(test$Z, score, target, t),
          stringsAsFactors = FALSE
        )
      )

      if (isTRUE(design$ld > 0)) {
        local_alt <- truth_alt_locus(alt, 5L)
        row$FP_loc <- sum(rejected & !local_alt)
        row$FDP_loc <- if (sum(rejected)) mean(!local_alt[rejected]) else 0
        for (window in c(1L, 2L, 5L)) {
          union_alt <- truth_alt(union_config(test$H, window), target, t)
          row[[paste0("FDP_union_", window)]] <- if (sum(rejected)) mean(!union_alt[rejected]) else 0
          row[[paste0("power_union_", window)]] <- if (sum(union_alt)) {
            sum(rejected & union_alt) / sum(union_alt)
          } else NA_real_
        }
      }
      if (length(extra)) row <- cbind(row, as.data.frame(extra, stringsAsFactors = FALSE))
      row
    }

    rows <- NULL
    safe_fit <- function(z, model, control) {
      tryCatch(csm_fit(z, model, starts = c("null", "shared"), control = control),
               error = function(e) list(converged = FALSE, message = conditionMessage(e)))
    }

    if (grepl("^ours", method_type)) {
      component_pattern <- if (grepl("_[AC]$", method_type)) {
        sub(".*_([AC])$", "\\1", method_type)
      } else "B"
      components <- switch(component_pattern,
                           A = c(1L, 1L),
                           B = c(2L, 1L),
                           C = c(2L, 2L))
      penalized <- grepl("_pen", method_type)
      lambda_mu <- 30
      lambda_weight <- 10
      label_suffix <- ""
      match <- regexec("_L([0-9]+)A([0-9]+)$", method_type)
      parsed <- regmatches(method_type, match)[[1L]]
      if (length(parsed)) {
        lambda_mu <- as.numeric(parsed[2L])
        lambda_weight <- as.numeric(parsed[3L])
        label_suffix <- sprintf("_L%gA%g", lambda_mu, lambda_weight)
      }
      control <- if (penalized) {
        csm_control(cores = cores, max_iter = max_iter,
                    penalty = lambda_mu, weight_penalty = lambda_weight)
      } else {
        csm_control(cores = cores, max_iter = max_iter, penalty = 0, weight_penalty = 0)
      }
      method_label <- paste0(if (penalized) "penalized" else "ML",
                             if (component_pattern == "B") "" else paste0("_", component_pattern),
                             label_suffix)

      if (grepl("_oracle$", method_type)) {
        null_rows <- rowSums(train$H) == 0
        oracle_cor <- stats::cor(train$Z[null_rows, , drop = FALSE])
        oracle_cor <- csmGmmCorr:::.project_cor(oracle_cor, 0.05)
        control <- csm_control(cores = cores, max_iter = max_iter,
                               penalty = lambda_mu, weight_penalty = lambda_weight,
                               fix_sigma = oracle_cor)
        method_label <- "penalized_nullSigma"
      }

      model <- csm_spec(K, M = csm_structure(K, components[1L], components[2L]), b_mu = 15)
      fit_time <- system.time(fit <- safe_fit(train$Z, model, control))[["elapsed"]]
      message(sprintf("replicate %d: %s fit %.1f min, converged=%s",
                      replicate, method_label, fit_time / 60, isTRUE(fit$converged)))
      fit_info <- list(
        converged = isTRUE(fit$converged),
        iterations = if (isTRUE(fit$converged)) sum(fit$runs$iterations) else NA_integer_,
        fit_min = fit_time / 60,
        gen_min = generated / 60,
        Sigma_minEig = if (isTRUE(fit$converged)) min(eigen(fit$params$Sigma, TRUE, TRUE)$values) else NA_real_,
        structure = component_pattern,
        loglik = if (isTRUE(fit$converged)) as.numeric(logLik(fit)) else NA_real_,
        df = if (isTRUE(fit$converged)) attr(logLik(fit), "df") else NA_integer_
      )

      for (tg in target) for (threshold in thresholds) {
        if (!isTRUE(fit$converged)) {
          rows <- rbind(rows, make_row(tg, threshold, method_label, rep(FALSE, N),
                                       c(fit_info, list(flagged = NA_integer_))))
          next
        }
        lfdr <- csm_lfdr(fit, test$Z, tg, threshold, cores = cores)
        selected <- csm_stepup(lfdr, q)
        extrapolated <- csm_extrapolation(fit, test$Z, tg, threshold,
                                          selected$rejected, selected$cutoff)
        rows <- rbind(
          rows,
          make_row(tg, threshold, method_label, selected$rejected,
                   c(fit_info, list(flagged = sum(extrapolated$flag, na.rm = TRUE))), lfdr),
          make_row(tg, threshold, paste0(method_label, "_adjusted"),
                   csm_adjust(test$Z, lfdr, tg, threshold, q)$rejected,
                   c(fit_info, list(flagged = NA_integer_)))
        )
      }
    } else if (grepl("_csm$", method_type)) {
      model <- csm_spec(K, M = csm_structure(K, 2L, 1L), b_mu = 15)
      plug_in <- csmGmmCorr:::.project_cor(stats::cor(train$Z), model$kappa)
      controls <- list(
        csmGmm = csm_control(cores = cores, max_iter = max_iter, penalty = 0,
                             weight_penalty = 0, fix_sigma = TRUE),
        csmGmm_cor = csm_control(cores = cores, max_iter = max_iter, penalty = 0,
                                 weight_penalty = 0, fix_sigma = plug_in),
        csmGmm_pen = csm_control(cores = cores, max_iter = max_iter, fix_sigma = TRUE)
      )
      for (method_label in names(controls)) {
        fit_time <- system.time(fit <- safe_fit(train$Z, model, controls[[method_label]]))[["elapsed"]]
        fit_info <- list(
          converged = isTRUE(fit$converged),
          fit_min = fit_time / 60,
          gen_min = generated / 60,
          structure = "B"
        )
        for (tg in target) for (threshold in thresholds) {
          if (!isTRUE(fit$converged)) {
            rows <- rbind(rows, make_row(tg, threshold, method_label, rep(FALSE, N), fit_info))
          } else {
            lfdr <- csm_lfdr(fit, test$Z, tg, threshold, cores = cores)
            selected <- csm_stepup(lfdr, q)
            rows <- rbind(rows, make_row(tg, threshold, method_label,
                                         selected$rejected, fit_info, lfdr))
          }
        }
      }
    } else {
      comparison <- switch(
        method_type,
        comp_pl = list(target = "pleiotropy", methods = list(
          QCH = function(z, t) run_qch_pleiotropy_R1(z, t, q, threads = cores),
          adaFilter = function(z, t) run_adafilter_pleiotropy_R1(z, t, q)
        )),
        comp_rp = list(target = "replication", methods = list(
          repfdr = function(z, t) run_repfdr_replication_R1(z, t, q),
          DirectionalPC = function(z, t) run_dirpc_replication_R1(z, t, q)
        )),
        comp_med = list(target = "mediation", methods = list(
          DEI = function(z, t) run_dei_mediation_R1(z, t, q),
          MediationTest = function(z, t) run_mediation_test_R1(z, t, q),
          PCPathSafeguard = function(z, t) run_mediation_safeguard_R1(z, t, q)
        )),
        stop("Unknown method type: ", method_type)
      )
      for (threshold in thresholds) for (method_label in names(comparison$methods)) {
        answer <- tryCatch(comparison$methods[[method_label]](test$Z, threshold),
                           error = function(e) list(error = conditionMessage(e)))
        if (!is.null(answer$error)) {
          rows <- rbind(rows, cbind(make_row(comparison$target, threshold, method_label,
                                             rep(FALSE, N)), error = answer$error))
        } else {
          rows <- rbind(rows, cbind(make_row(comparison$target, threshold, method_label,
                                             answer$rejected, score = answer$score), error = ""))
        }
      }
    }

    rows$total_min <- (proc.time()[["elapsed"]] - started) / 60
    utils::write.csv(rows, destination, row.names = FALSE)
    message("done: ", destination)
    invisible(NULL)
  }

  replicates <- seq.int(task$rep_from, task$rep_to)
  maximum <- as.integer(Sys.getenv("SIM_MAX_REPS", "0"))
  if (maximum > 0L) replicates <- utils::head(replicates, maximum)
  for (replicate in replicates) evaluate_one(replicate)
  invisible(NULL)
}
