RD_ANALYSES <- list(
  glgc_replication = list(
    target = "replication", prep = "glgc_replication",
    cols = c("EUR_female", "EUR_male", "EAS_female", "EAS_male", "AFR_female", "AFR_male"),
    holdout = c("HIS_female", "HIS_male", "SAS_female", "SAS_male")
  ),
  ukb_pleiotropy = list(
    target = "pleiotropy", prep = "ukb_pleiotropy",
    cols = c("BMI", "SBP", "TG", "HDL", "HbA1c")
  ),
  glgc_pleiotropy = list(
    target = "pleiotropy", prep = "glgc_pleiotropy",
    cols = c("HDL", "LDL", "TC", "logTG", "nonHDL")
  )
)

RD_JOBS <- data.frame(
  analysis = c(rep("glgc_replication", 5L), rep("ukb_pleiotropy", 4L), rep("glgc_pleiotropy", 4L)),
  variant = c("main", "indep", "robust", "main", "indep", "main", "indep", "robust", "main", "main", "indep", "robust", "main"),
  a_mu = c(0.5, 0.5, 0.5, 1.5, 1.5, 0.5, 0.5, 0.5, 1.5, 0.5, 0.5, 0.5, 1.5),
  stringsAsFactors = FALSE
)

rd_label <- function(variant, a_mu) paste0(variant, if (a_mu != 0.5) sprintf("_amu%g", 10 * a_mu) else "")

rd_read_gz <- function(path, ...) {
  if (grepl("\\.(gz|bgz)$", path)) {
    data.table::fread(cmd = paste("gzip -cd", shQuote(path)), ...)
  } else data.table::fread(path, ...)
}

rd_cap_z <- function(z, cap = 12) {
  storage.mode(z) <- "double"
  pmin(pmax(z, -cap), cap)
}

rd_load <- function(analysis, cap = 12, smoke = FALSE, extra = NULL,
                    root = getOption("csmgmm.reproduce.root", getwd())) {
  spec <- RD_ANALYSES[[analysis]]
  if (is.null(spec)) stop("Unknown real-data analysis: ", analysis)
  prepared <- file.path(root, "Data", "prepared", spec$prep)
  manifest_path <- file.path(prepared, "chromosome_manifest.tsv")
  if (!file.exists(manifest_path)) stop("Prepared data not found: ", prepared)
  fit_sample <- rd_read_gz(file.path(prepared, "fit_sample.tsv.gz"))
  if (smoke) {
    set.seed(1)
    fit_sample <- fit_sample[sample.int(nrow(fit_sample), min(nrow(fit_sample), 20000L))]
  }
  statistic_cols <- spec$cols
  Zfit <- rd_cap_z(as.matrix(fit_sample[, ..statistic_cols]), cap)
  manifest <- data.table::fread(manifest_path)
  if (smoke) manifest <- manifest[as.character(chr) %in% "22"]
  parts <- lapply(seq_len(nrow(manifest)), function(i) {
    d <- rd_read_gz(file.path(prepared, manifest$path[i]),
                    select = unique(c("variant_id", "chr", "pos_b37", spec$cols, extra)))
    if (nrow(d) == manifest$n_rows[i]) d else stop("Row-count mismatch in ", manifest$path[i])
  })
  all <- data.table::rbindlist(parts, fill = TRUE)
  keep <- !all$variant_id %in% fit_sample$variant_id
  all <- all[keep]
  Z <- rd_cap_z(as.matrix(all[, ..statistic_cols]), cap)
  info_cols <- intersect(c("variant_id", "chr", "pos_b37", extra), names(all))
  list(Zfit = Zfit, Z = Z, info = all[, ..info_cols],
       n_fit_excluded = sum(!keep), columns = spec$cols, target = spec$target)
}

robust_cor <- function(Z) {
  scale <- apply(Z, 2L, function(x) stats::IQR(x) / 1.349)
  K <- ncol(Z)
  out <- diag(K)
  for (i in seq_len(K - 1L)) for (j in (i + 1L):K) {
    a <- Z[, i] / scale[i]
    b <- Z[, j] / scale[j]
    out[i, j] <- out[j, i] <- (stats::var(a + b) - stats::var(a - b)) /
      (stats::var(a + b) + stats::var(a - b))
  }
  dimnames(out) <- list(colnames(Z), colnames(Z))
  out
}

project_cor <- function(R, kappa = 0.05, margin = 1e-3) {
  R <- (R + t(R)) / 2
  for (iteration in seq_len(1000L)) {
    eig <- eigen(R, symmetric = TRUE)
    if (min(eig$values) >= kappa + margin / 2) break
    R <- eig$vectors %*% diag(pmax(eig$values, kappa + margin)) %*% t(eig$vectors)
    R <- stats::cov2cor(R)
  }
  R <- (R + t(R)) / 2
  diag(R) <- 1
  R
}

reference_cors <- function(Z) {
  keep <- rowSums(abs(Z) < 2) == ncol(Z)
  list(
    pearson = stats::cor(Z),
    pearson_abs_lt2 = stats::cor(Z[keep, , drop = FALSE]),
    robust = robust_cor(Z),
    iqr_sd = apply(Z, 2L, function(v) stats::IQR(v) / 1.349)
  )
}

rd_majorant <- function(Z, score, target, t) {
  if (is.loaded("rd_majorant")) {
    return(.Call("rd_majorant", Z, as.double(score), target,
                 as.integer(t), 32L))
  }
  csmGmmCorr::csm_adjust(Z, score, target, t, all = TRUE)$adjusted
}
rd_incon_dec <- function(Z, rejected, target, t) {
  if (is.loaded("rd_incon_dec")) {
    violating <- .Call("rd_incon_dec", Z, as.logical(rejected), target,
                       as.integer(t), 32L)
    return(list(violations = sum(violating), violating = violating,
                exact = TRUE))
  }
  score <- as.numeric(!as.logical(rejected))
  diagnostic <- csmGmmCorr::csm_incongruity(
    Z, score, target, t, max_n = nrow(Z)
  )
  list(violations = diagnostic$violations, violating = NULL, exact = TRUE)
}

clump_leads <- function(chr, pos, score, selected, window = 5e5) {
  idx <- which(selected)
  if (!length(idx)) return(integer())
  idx <- idx[order(score[idx], pos[idx])]
  leaders <- integer()
  occupied <- list()
  for (i in idx) {
    chromosome <- as.character(chr[i])
    positions <- occupied[[chromosome]]
    if (length(positions) && any(abs(positions - pos[i]) <= window)) next
    occupied[[chromosome]] <- sort(c(positions, pos[i]))
    leaders <- c(leaders, i)
  }
  leaders
}

clump_count <- function(chr, pos, score, selected, window = 5e5) {
  length(clump_leads(chr, pos, score, selected, window))
}

run_realdata <- function(analysis, variant = "main", a_mu = 0.5,
                         out_root = file.path(getOption("csmgmm.reproduce.root", getwd()), "Fig4", "output"),
                         cores = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "8")),
                         q = 0.1, cap = 12, ts = 1:2, smoke = identical(Sys.getenv("SMOKE"), "1")) {
  if (!analysis %in% names(RD_ANALYSES)) stop("Unknown analysis")
  if (!variant %in% c("main", "indep", "robust")) stop("Unknown covariance variant.")
  setup <- RD_ANALYSES[[analysis]]
  label <- rd_label(variant, a_mu)
  out <- file.path(out_root, paste0(analysis, "_", label))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(out, "DONE")) && Sys.getenv("RERUN", "0") != "1") return(invisible(out))
  dat <- rd_load(analysis, cap, smoke)
  ref <- reference_cors(dat$Zfit)
  fixed <- switch(variant, main = FALSE, indep = TRUE, robust = project_cor(ref$robust), FALSE)
  M <- csmGmmCorr::csm_structure(ncol(dat$Z), 2L, 1L)
  spec <- csmGmmCorr::csm_spec(ncol(dat$Z), M, S = ncol(dat$Z), a_mu = a_mu, b_mu = 15)
  control <- csmGmmCorr::csm_control(max_iter = if (smoke) 300L else 5000L,
                                    penalty = 30, weight_penalty = 10,
                                    fix_sigma = fixed, cores = cores)
  start <- proc.time()[3L]
  fit <- csmGmmCorr::csm_fit(dat$Zfit, spec, starts = c("null", "shared"), control = control)
  fit_minutes <- (proc.time()[3L] - start) / 60
  saveRDS(fit, file.path(out, "fit.rds"))
  fit_runs <- data.table::as.data.table(fit$runs)
  fit_runs[, `:=`(analysis = analysis, variant = label)]
  data.table::setcolorder(fit_runs, c("analysis", "variant", setdiff(names(fit_runs), c("analysis", "variant"))))
  data.table::fwrite(fit_runs, file.path(out, "fit_runs.tsv"), sep = "\t")

  matrix_rows <- function(matrix, kind) {
    index <- which(upper.tri(matrix), arr.ind = TRUE)
    data.table::data.table(
      analysis = analysis, variant = label, kind = kind,
      stat1 = dat$columns[index[, 1L]], stat2 = dat$columns[index[, 2L]],
      value = matrix[index]
    )
  }
  fitted_sigma <- if (!is.null(fit$params)) fit$params$Sigma else matrix(NA_real_, ncol(dat$Z), ncol(dat$Z))
  fixed_input <- if (is.matrix(fixed)) fixed else if (isTRUE(fixed)) diag(ncol(dat$Z)) else NULL
  sigma_parts <- list(
    matrix_rows(fitted_sigma, "fitted_Sigma"),
    matrix_rows(ref$pearson, "pearson"),
    matrix_rows(ref$pearson_abs_lt2, "pearson_abs_lt2"),
    matrix_rows(ref$robust, "robust")
  )
  if (!is.null(fixed_input)) sigma_parts[[length(sigma_parts) + 1L]] <- matrix_rows(fixed_input, "fixed_Sigma_input")
  data.table::fwrite(data.table::rbindlist(sigma_parts), file.path(out, "sigma.tsv"), sep = "\t")
  data.table::fwrite(
    data.table::data.table(analysis = analysis, variant = label,
                           stat = dat$columns, iqr_sd = ref$iqr_sd),
    file.path(out, "iqr_sd.tsv"), sep = "\t"
  )
  if (!is.null(fit$params)) {
    mu <- fit$params$mu
    colnames(mu) <- paste0("mu_", dat$columns)
    support <- data.table::data.table(
      analysis = analysis, variant = label,
      support = apply(fit$spec$grp_A, 1L, paste0, collapse = ""),
      size = fit$spec$grp_size, component = fit$spec$grp_m,
      weight = fit$params$varpi
    )
    data.table::fwrite(cbind(support, data.table::as.data.table(mu)),
                       file.path(out, "support_weights.tsv"), sep = "\t")
  }

  set.seed(20261006)
  sample_index <- sort(sample.int(nrow(dat$Z), min(nrow(dat$Z), 20000L)))
  plot_sample <- data.table::data.table(dat$info[sample_index], dat$Z[sample_index, , drop = FALSE])
  results <- list()
  for (t in ts) {
    if (!isTRUE(fit$converged)) break
    infer_start <- proc.time()[3L]
    lfdr <- csmGmmCorr::csm_lfdr(fit, dat$Z, setup$target, t, cores = cores)
    step <- csmGmmCorr::csm_stepup(lfdr, q)
    adjusted <- rd_majorant(dat$Z, lfdr, setup$target, t)
    if (any(adjusted + 1e-12 < lfdr)) stop("Majorant score fell below raw lfdr.")
    adjusted_step <- csmGmmCorr::csm_stepup(adjusted, q)
    decision_raw <- rd_incon_dec(dat$Z, step$rejected, setup$target, t)
    decision_adjusted <- rd_incon_dec(dat$Z, adjusted_step$rejected, setup$target, t)
    incon_raw <- csmGmmCorr::csm_incongruity(dat$Z, lfdr, setup$target, t, max_n = 5000L)
    incon_adjusted <- csmGmmCorr::csm_incongruity(dat$Z, adjusted, setup$target, t, max_n = 5000L)
    extrapolation <- if (step$n_rejected) {
      csmGmmCorr::csm_extrapolation(
        fit, dat$Z, setup$target, t,
        rejected = step$rejected, cutoff = step$cutoff
      )$flag
    } else rep(NA, nrow(dat$Z))
    selected_iterations <- fit$runs$iterations[fit$runs$start == fit$selected]
    results[[length(results) + 1L]] <- data.frame(
      analysis, variant = label, target = setup$target, K = ncol(dat$Z), t, q,
      n_fit = nrow(dat$Zfit), n_test = nrow(dat$Z), converged = fit$converged,
      selected_start = fit$selected,
      iterations = if (length(selected_iterations)) selected_iterations[1L] else NA_integer_,
      loglik = fit$loglik, fit_minutes = round(fit_minutes, 2),
      R_raw = step$n_rejected, R_adjusted = adjusted_step$n_rejected,
      R_both = sum(step$rejected & adjusted_step$rejected),
      frac_raw = step$n_rejected / nrow(dat$Z), cutoff_raw = step$cutoff,
      avg_lfdr_raw = step$average,
      loci_raw = clump_count(dat$info$chr, dat$info$pos_b37, lfdr, step$rejected),
      loci_adjusted = clump_count(dat$info$chr, dat$info$pos_b37, adjusted, adjusted_step$rejected),
      n_adjusted_up = sum(adjusted > lfdr + 1e-12),
      incon_score_raw = incon_raw$violations,
      incon_score_adjusted = incon_adjusted$violations,
      incon_decision_raw = decision_raw$violations,
      incon_decision_adjusted = decision_adjusted$violations,
      extrapolation_flags = sum(extrapolation, na.rm = TRUE),
      inference_minutes = round((proc.time()[3L] - infer_start) / 60, 2)
    )
    violating <- rep(FALSE, nrow(dat$Z))
    if (!is.null(decision_raw$violating)) violating <- violating | decision_raw$violating
    if (!is.null(decision_adjusted$violating)) violating <- violating | decision_adjusted$violating
    keep <- step$rejected | adjusted_step$rejected | violating
    discoveries <- data.table::data.table(
      dat$info[keep], dat$Z[keep, , drop = FALSE],
      lfdr = lfdr[keep], lfdr_adjusted = adjusted[keep],
      reject_raw = step$rejected[keep], reject_adjusted = adjusted_step$rejected[keep],
      incongruous_decision = violating[keep], extrapolation = extrapolation[keep]
    )
    data.table::fwrite(discoveries, file.path(out, sprintf("discoveries_t%d.tsv.gz", t)), sep = "\t")
    plot_sample[[paste0("lfdr_t", t)]] <- lfdr[sample_index]
    plot_sample[[paste0("lfdr_adjusted_t", t)]] <- adjusted[sample_index]
  }
  data.table::fwrite(plot_sample, file.path(out, "plot_sample.tsv.gz"), sep = "\t")
  if (length(results)) data.table::fwrite(data.table::rbindlist(results), file.path(out, "summary.tsv"), sep = "\t")
  writeLines("DONE", file.path(out, "DONE"))
  invisible(out)
}

rd_fit_dirs <- function(out_root) {
  file.path(out_root, paste0(RD_JOBS$analysis, "_", mapply(rd_label, RD_JOBS$variant, RD_JOBS$a_mu)))
}

rd_collect <- function(out_root = file.path(getOption("csmgmm.reproduce.root", getwd()), "Fig4", "output")) {
  dirs <- rd_fit_dirs(out_root)
  read_named <- function(file) {
    rows <- lapply(dirs[file.exists(file.path(dirs, file))], function(d) {
      x <- data.table::fread(file.path(d, file))
      x[, fit_label := basename(d)]
      x
    })
    data.table::rbindlist(rows, fill = TRUE)
  }
  list(summary = read_named("summary.tsv"), sigma = read_named("sigma.tsv"),
       fit_runs = read_named("fit_runs.tsv"), iqr_sd = read_named("iqr_sd.tsv"))
}

MEDLDL_COLUMNS <- c("alpha1", "beta1", "alpha2", "beta2", "alpha3", "beta3")
MEDLDL_JOBS <- rbind(
  expand.grid(set = "noshare", variant = c("main", "indep", "robust"), a_mu = c(.5, 1.5), stringsAsFactors = FALSE),
  expand.grid(set = "nearest", variant = "main", a_mu = c(.5, 1.5), stringsAsFactors = FALSE),
  expand.grid(set = "noshare_dedup", variant = c("main", "indep", "robust"), a_mu = c(.5, 1.5), stringsAsFactors = FALSE),
  expand.grid(set = "noshare_cond", variant = c("main", "indep", "robust"), a_mu = c(.5, 1.5), stringsAsFactors = FALSE),
  expand.grid(set = "noshare", variant = c("main", "indep", "robust"), a_mu = 3, stringsAsFactors = FALSE)
)

medldl_dedup <- function(d) {
  gene_cols <- paste0("gene", 1:3)
  key <- apply(d[, ..gene_cols], 1L, function(x) paste(sort(x), collapse = ","))
  copy <- data.table::copy(d)
  copy[, triple := key]
  data.table::setorder(copy, triple, chr, pos_b38)
  copy[, .SD[1L], by = triple][, triple := NULL]
}

run_mediation_ldl <- function(set, variant = "main", a_mu = 0.5,
                              out_root = file.path(getOption("csmgmm.reproduce.root", getwd()), "Fig4", "output", "mediation_ldl"),
                              cores = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "4")),
                              q = 0.1, cap = 12, ts = 1:2) {
  if (!set %in% c("noshare", "nearest", "noshare_dedup", "noshare_cond")) stop("Unknown mediation set.")
  if (!variant %in% c("main", "indep", "robust")) stop("Unknown covariance variant.")
  input_set <- sub("_dedup$", "", set)
  path <- file.path(getOption("csmgmm.reproduce.root", getwd()), "Data", "mediation",
                    sprintf("mediation_ldl_liver_%s.tsv.gz", input_set))
  d <- rd_read_gz(path)
  if (grepl("_dedup$", set)) d <- medldl_dedup(d)
  gene_cols <- paste0("gene", 1:3)
  required <- c(MEDLDL_COLUMNS, "chr", "pos_b38", gene_cols)
  missing <- setdiff(required, names(d))
  if (length(missing)) stop("Mediation input is missing: ", paste(missing, collapse = ", "))
  d <- d[stats::complete.cases(d[, ..required])]
  Z <- rd_cap_z(as.matrix(d[, ..MEDLDL_COLUMNS]), cap)
  label <- paste0(set, "_", rd_label(variant, a_mu))
  out <- file.path(out_root, label)
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  if (file.exists(file.path(out, "DONE")) && Sys.getenv("RERUN", "0") != "1") {
    return(invisible(out))
  }
  ref <- reference_cors(Z)
  fixed <- switch(variant, main = FALSE, indep = TRUE, robust = project_cor(ref$robust), FALSE)
  M <- csmGmmCorr::csm_structure(ncol(Z), 2L, 1L)
  spec <- csmGmmCorr::csm_spec(ncol(Z), M, S = ncol(Z), a_mu = a_mu, b_mu = 15)
  control <- csmGmmCorr::csm_control(max_iter = 5000L, penalty = 30, weight_penalty = 10,
                                    fix_sigma = fixed, cores = cores)
  started <- proc.time()[3L]
  fit <- csmGmmCorr::csm_fit(Z, spec, starts = c("null", "shared"), control = control)
  fit_minutes <- (proc.time()[3L] - started) / 60
  saveRDS(fit, file.path(out, "fit.rds"))

  matrix_rows <- function(matrix, kind) {
    index <- which(upper.tri(matrix), arr.ind = TRUE)
    data.table::data.table(
      set = set, variant = rd_label(variant, a_mu), kind = kind,
      stat1 = colnames(Z)[index[, 1L]], stat2 = colnames(Z)[index[, 2L]],
      value = matrix[index]
    )
  }
  sigma <- if (!is.null(fit$params)) fit$params$Sigma else matrix(NA_real_, ncol(Z), ncol(Z))
  sigma_rows <- data.table::rbindlist(list(
    matrix_rows(sigma, "fitted"),
    matrix_rows(ref$pearson, "pearson"),
    matrix_rows(ref$pearson_abs_lt2, "pearson_abs_lt2"),
    matrix_rows(ref$robust, "robust")
  ))
  data.table::fwrite(sigma_rows, file.path(out, "sigma.tsv"), sep = "\t")
  data.table::fwrite(
    data.table::data.table(set = set, variant = rd_label(variant, a_mu),
                           stat = colnames(Z), iqr_sd = ref$iqr_sd),
    file.path(out, "iqr_sd.tsv"), sep = "\t"
  )
  if (!is.null(fit$params)) {
    mu <- fit$params$mu
    colnames(mu) <- paste0("mu_", colnames(Z))
    support <- data.table::data.table(
      set = set, variant = rd_label(variant, a_mu),
      support = apply(fit$spec$grp_A, 1L, paste0, collapse = ""),
      size = fit$spec$grp_size, component = fit$spec$grp_m,
      weight = fit$params$varpi
    )
    support <- cbind(support, data.table::as.data.table(mu))
    data.table::fwrite(support, file.path(out, "support_weights.tsv"), sep = "\t")
  }

  triple_count <- function(selected) {
    if (!any(selected)) return(0L)
    length(unique(apply(d[selected, ..gene_cols], 1L, function(x) paste(sort(x), collapse = ":"))))
  }
  rows <- list()
  discoveries <- list()
  for (t in ts) {
    if (!isTRUE(fit$converged)) break
    lfdr <- csmGmmCorr::csm_lfdr(fit, Z, "mediation", t, cores = cores)
    raw <- csmGmmCorr::csm_stepup(lfdr, q)
    adjusted <- rd_majorant(Z, lfdr, "mediation", t)
    adj <- csmGmmCorr::csm_stepup(adjusted, q)
    extrapolation <- if (raw$n_rejected) {
      csmGmmCorr::csm_extrapolation(
        fit, Z, "mediation", t, rejected = raw$rejected,
        cutoff = raw$cutoff
      )$flag
    } else rep(NA, nrow(Z))
    rows[[length(rows) + 1L]] <- data.frame(
      set, variant = rd_label(variant, a_mu), a_mu, t, n_features = nrow(Z),
      n_triples = nrow(medldl_dedup(d)), converged = fit$converged,
      loglik = fit$loglik, fit_minutes = round(fit_minutes, 2),
      null_weight = sum(fit$params$varpi[fit$spec$grp_size == 0L]),
      R_raw = raw$n_rejected, R_adjusted = adj$n_rejected,
      triples_raw = triple_count(raw$rejected),
      loci_raw = clump_count(as.character(d$chr), d$pos_b38, lfdr, raw$rejected),
      loci_adjusted = clump_count(as.character(d$chr), d$pos_b38, lfdr, adj$rejected),
      extrapolation_flags = sum(extrapolation, na.rm = TRUE),
      incon_raw = csmGmmCorr::csm_incongruity(Z, lfdr, "mediation", t)$violations,
      incon_adjusted = csmGmmCorr::csm_incongruity(Z, adjusted, "mediation", t)$violations
    )
    keep <- raw$rejected | adj$rejected
    discoveries[[length(discoveries) + 1L]] <- data.table::data.table(
      t = t, d[keep], Z[keep, , drop = FALSE], lfdr = lfdr[keep], lfdr_adjusted = adjusted[keep],
      reject_raw = raw$rejected[keep], reject_adjusted = adj$rejected[keep]
    )
  }
  if (length(rows)) data.table::fwrite(data.table::rbindlist(rows), file.path(out, "summary.tsv"), sep = "\t")
  if (length(discoveries)) data.table::fwrite(data.table::rbindlist(discoveries, fill = TRUE), file.path(out, "discoveries.tsv.gz"), sep = "\t")
  writeLines("DONE", file.path(out, "DONE"))
  invisible(out)
}

medldl_fit_dirs <- function(out_root = file.path(getOption("csmgmm.reproduce.root", getwd()),
                                                  "Fig4", "output", "mediation_ldl")) {
  file.path(out_root, paste0(MEDLDL_JOBS$set, "_",
                             mapply(rd_label, MEDLDL_JOBS$variant, MEDLDL_JOBS$a_mu)))
}

run_mediation_all_t <- function(index, out_root = file.path(getOption("csmgmm.reproduce.root", getwd()),
                                                             "Fig4", "output", "mediation_ldl"),
                                data_root = file.path(getOption("csmgmm.reproduce.root", getwd()),
                                                      "Data", "mediation"),
                                cores = as.integer(Sys.getenv("LSB_DJOB_NUMPROC", "4")),
                                q = 0.1, cap = 12) {
  dirs <- medldl_fit_dirs(out_root)
  available <- dirs[file.exists(file.path(dirs, "fit.rds"))]
  available <- sort(available)
  index <- as.integer(index)
  if (is.na(index) || index < 1L || index > length(available)) return(invisible(NULL))
  out <- available[index]
  label <- basename(out)
  matched <- which(basename(dirs) == label)
  if (length(matched) != 1L) stop("Cannot identify mediation fit: ", label)
  job <- MEDLDL_JOBS[matched, ]
  input_set <- sub("_dedup$", "", job$set)
  d <- rd_read_gz(file.path(data_root, sprintf("mediation_ldl_liver_%s.tsv.gz", input_set)))
  if (grepl("_dedup$", job$set)) d <- medldl_dedup(d)
  required <- c(MEDLDL_COLUMNS, "chr", "pos_b38", paste0("gene", 1:3))
  d <- d[stats::complete.cases(d[, ..required])]
  Z <- rd_cap_z(as.matrix(d[, ..MEDLDL_COLUMNS]), cap)
  fit <- readRDS(file.path(out, "fit.rds"))
  if (!isTRUE(fit$converged)) stop("Saved mediation fit did not converge: ", label)

  magnitude <- abs(Z)
  strong_path <- (magnitude[, 1L] >= 4 & magnitude[, 2L] >= 4) |
    (magnitude[, 3L] >= 4 & magnitude[, 4L] >= 4) |
    (magnitude[, 5L] >= 4 & magnitude[, 6L] >= 4)
  n_strong <- rowSums(magnitude >= 4)
  max_alpha <- apply(magnitude[, c(1L, 3L, 5L), drop = FALSE], 1L, max)
  rows <- lapply(1:5, function(t) {
    lfdr <- csmGmmCorr::csm_lfdr(fit, Z, "mediation", t, cores = cores)
    raw <- csmGmmCorr::csm_stepup(lfdr, q)
    adjusted <- rd_majorant(Z, lfdr, "mediation", t)
    adj <- csmGmmCorr::csm_stepup(adjusted, q)
    leaders <- clump_leads(as.character(d$chr), d$pos_b38, lfdr,
                           rep(TRUE, nrow(d)))
    locus <- integer(nrow(d))
    for (chromosome in unique(as.character(d$chr))) {
      members <- which(as.character(d$chr) == chromosome)
      chromosome_leads <- leaders[as.character(d$chr[leaders]) == chromosome]
      if (!length(chromosome_leads)) next
      locus[members] <- vapply(d$pos_b38[members], function(position) {
        chromosome_leads[which.min(abs(d$pos_b38[chromosome_leads] - position))]
      }, integer(1))
    }
    data.table::data.table(
      fit = label, t = t, n = nrow(Z), R_raw = raw$n_rejected,
      R_adjusted = adj$n_rejected,
      strong_target_data = sum(strong_path),
      strong_target_loci_data = if (any(strong_path))
        clump_count(as.character(d$chr), d$pos_b38, -n_strong, strong_path) else 0L,
      strong_target_rejections = sum(raw$rejected & strong_path),
      all_alpha_lt2 = sum(raw$rejected & max_alpha < 2),
      loci = length(unique(locus[raw$rejected])),
      loci_with_strong = length(unique(locus[raw$rejected & strong_path]))
    )
  })
  data.table::fwrite(data.table::rbindlist(rows), file.path(out, "all_t.tsv"), sep = "\t")
  invisible(file.path(out, "all_t.tsv"))
}
