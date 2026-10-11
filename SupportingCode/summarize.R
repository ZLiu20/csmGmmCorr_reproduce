read_sim <- function(folder, study, refresh = FALSE, root = here::here()) {
  output_dir <- file.path(root, folder, "output")
  cache <- file.path(output_dir, sprintf("all_results_%s.csv.gz", study))
  if (file.exists(cache) && !refresh) return(data.table::fread(cache))
  files <- list.files(file.path(output_dir, study), pattern = "[.]csv$",
                      recursive = TRUE, full.names = TRUE)
  if (!length(files)) stop("No simulation results in ", file.path(folder, "output", study))
  message(sprintf("Reading %d files for %s/%s", length(files), folder, study))
  result <- data.table::rbindlist(lapply(files, data.table::fread), fill = TRUE)
  data.table::fwrite(result, cache)
  result
}

sim_progress <- function(folder = NULL, root = here::here()) {
  arrays <- sim_arrays()
  if (!is.null(folder)) arrays <- arrays[vapply(arrays, function(x) x$folder %in% folder, logical(1))]
  do.call(rbind, lapply(arrays, function(x) {
    tasks <- sim_task_table(x$scen, x$designs, x$reps)
    output <- file.path(root, x$folder, "output", x$study, x$type)
    data.frame(array = x$name, folder = x$folder, tasks = nrow(tasks),
               expected = sum(tasks$rep_to - tasks$rep_from + 1L),
               found = length(list.files(output, pattern = "[.]csv$")))
  }))
}

design_means <- function(d, by = c("target", "K", "t", "method", "design", "axis", "x")) {
  d <- data.table::as.data.table(d)
  d[, .(FDP = mean(FDP), FDP_se = stats::sd(FDP) / sqrt(.N), n = .N,
        power = mean(power, na.rm = TRUE), R = mean(R),
        incon_dec = mean(incon_dec, na.rm = TRUE),
        incon_score = mean(incon_score, na.rm = TRUE)), by = by]
}

overall_means <- function(d, by = c("target", "K", "t", "method")) {
  d <- data.table::as.data.table(d)
  d[, .(FDP = mean(FDP), FDPmax = max(FDP), power = mean(power, na.rm = TRUE),
        ndes = .N, incon_dec = mean(incon_dec, na.rm = TRUE),
        incon_score = mean(incon_score, na.rm = TRUE)), by = by]
}

three_designs <- function(d) {
  d <- data.table::as.data.table(d)
  d[design %in% c("effect_01", "effect_05") |
      (scen != "med" & design == "tau1_002") |
      (scen == "med" & design == "tau1_004")]
}

f3 <- function(x, digits = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
fcomma <- function(x) ifelse(is.na(x), "--", gsub(",", "{,}", formatC(x, format = "d", big.mark = ",")))
tex_escape <- function(x) gsub("_", "\\\\_", x)
write_tex <- function(lines, file) {
  writeLines(lines, file)
  message("wrote ", file)
}

SOURCE_ROWS <- list(
  c("pen", "Correlated csmGmm (penalized)"),
  c(NA, "\\textit{Estimation}"),
  c("ML", "Correlated csmGmm, maximum likelihood"),
  c("csm", "csmGmm ($\\Sigma=I$), maximum likelihood"),
  c("csm_cor", "csmGmm, $\\Sigma$ fixed at sample correlation"),
  c("csm_pen", "csmGmm ($\\Sigma=I$), penalized"),
  c("split", "Correlated csmGmm, independent test set"),
  c(NA, "\\textit{Data-generating process}"),
  c("sym", "Symmetric statistics (intercept 0)"),
  c("quant", "Quantitative outcomes"),
  c("exch", "Exchangeable correlation 0.4"),
  c("ld", "Linkage disequilibrium (variant level)"),
  c(NA, "\\textit{Strong correlation ($\\rho=0.8$, fifth effect-size design)}"),
  c("rho08", "Correlated csmGmm"),
  c("rho08_oracle", "Correlated csmGmm, $\\Sigma$ fixed at null correlation"),
  c("rho08_csm", "csmGmm ($\\Sigma=I$), maximum likelihood")
)

sources_means <- function(value = "FDP", root = here::here()) {
  average <- function(d, method) {
    if (is.null(d)) return(NULL)
    d <- data.table::as.data.table(d)[method == ..method]
    if (!nrow(d)) return(NULL)
    g <- d[, .(v = mean(get(value), na.rm = TRUE)), by = .(target, K, t, design)]
    g[, .(v = mean(v)), by = .(target, K, t)]
  }
  optional <- function(folder, study) tryCatch(read_sim(folder, study, root = root),
                                                error = function(e) { message(conditionMessage(e)); NULL })
  main <- read_sim("Fig2", "main", root = root)
  sym <- read_sim("Fig2", "sym", root = root)
  alternatives <- optional("Tab2", "main")
  rho <- optional("Fig3", "rho")
  oracle <- optional("Tab2", "rho")
  rho08 <- function(d) if (is.null(d)) NULL else d[design == "rho_080"]
  td <- three_designs
  pieces <- list(
    pen = average(td(main), "penalized"),
    ML = average(td(main), "ML"),
    csm = average(td(alternatives), "csmGmm"),
    csm_cor = average(td(alternatives), "csmGmm_cor"),
    csm_pen = average(td(alternatives), "csmGmm_pen"),
    split = average(td(optional("Tab2", "split")), "penalized"),
    sym = average(td(sym), "penalized"),
    quant = average(td(optional("Tab2", "quant")), "penalized"),
    exch = average(td(optional("Tab2", "exch")), "penalized"),
    ld = average(td(optional("Tab2", "ld")), "penalized"),
    rho08 = average(rho08(rho), "penalized"),
    rho08_oracle = average(rho08(oracle), "penalized_nullSigma"),
    rho08_csm = average(rho08(rho), "csmGmm")
  )
  answer <- unique(main[, .(target, K, t)])
  for (name in names(pieces)) {
    part <- pieces[[name]]
    if (is.null(part)) answer[, (name) := NA_real_]
    else answer <- merge(answer, setnames(part, "v", name), by = c("target", "K", "t"), all.x = TRUE)
  }
  data.table::setorder(answer, target, K, t)
  answer
}

sources_tex <- function(table, threshold, label, caption) {
  lines <- c(
    "\\begin{table}[!htbp]", "\\centering",
    sprintf("\\caption{%s}\\label{%s}", caption, label), "\\small",
    "\\begin{tabular}{l cc cc cc}", "\\toprule",
    " & \\multicolumn{2}{c}{Pleiotropy} & \\multicolumn{2}{c}{Replication} & \\multicolumn{2}{c}{Mediation} \\\\",
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
    " & $K=4$ & $K=6$ & $K=4$ & $K=6$ & $K=4$ & $K=6$ \\\\ \\midrule"
  )
  for (entry in SOURCE_ROWS) {
    key <- entry[1L]
    if (is.na(key)) {
      lines <- c(lines, sprintf("\\multicolumn{7}{l}{%s}\\\\", entry[2L]))
      next
    }
    cells <- unlist(lapply(c("pleiotropy", "replication", "mediation"), function(target) {
      sapply(c(4L, 6L), function(K) {
        value <- table[target == ..target & K == ..K & t == threshold][[key]]
        f3(if (length(value)) value else NA_real_)
      })
    }))
    lines <- c(lines, paste(c(entry[2L], cells), collapse = " & "), "\\\\")
  }
  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

metric_tex <- function(d, metric, caption, label,
                       group = c("target", "K", "t", "method")) {
  d <- data.table::as.data.table(d)
  stopifnot(metric %in% names(d))
  z <- d[, .(mean = mean(get(metric), na.rm = TRUE),
             se = stats::sd(get(metric), na.rm = TRUE) / sqrt(.N), n = .N), by = group]
  head <- paste(group, collapse = " & ")
  align <- paste0("l", paste(rep("c", length(group) - 1L + 3L), collapse = ""))
  lines <- c("\\begin{table}[!htbp]", "\\centering",
             sprintf("\\caption{%s}\\label{%s}", caption, label), "\\small",
             sprintf("\\begin{tabular}{%s}", align), "\\toprule",
             sprintf("%s & Mean & MCSE & Replicates \\\\ \\midrule", head))
  for (i in seq_len(nrow(z))) {
    cells <- c(vapply(z[i, ..group], as.character, character(1L)),
               f3(z$mean[i]), f3(z$se[i]), as.character(z$n[i]))
    lines <- c(lines, paste(cells, collapse = " & "), "\\\\")
  }
  c(lines, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}

write_metric_table <- function(d, metric, file, caption, label,
                               group = c("target", "K", "t", "method")) {
  write_tex(metric_tex(d, metric, caption, label, group), file)
  invisible(file)
}
