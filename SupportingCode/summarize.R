# Reading and summarizing the simulation output (used by the plot_* and make_* scripts).

# All replicate files of one study in one folder: <folder>/output/<study>/<type>/*.csv.
# The combined table is cached in <folder>/output/all_results_<study>.csv.gz; set
# refresh = TRUE after new replicates have been written.
read_sim <- function(folder, study, refresh = FALSE, root = here::here()) {
  od <- file.path(root, folder, "output")
  cache <- file.path(od, sprintf("all_results_%s.csv.gz", study))
  if (file.exists(cache) && !refresh) return(data.table::fread(cache))
  f <- list.files(file.path(od, study), pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
  if (!length(f)) stop("no results in ", file.path(od, study), "; run the arrays of ", folder, " first")
  message(sprintf("reading %d files of study %s in %s", length(f), study, folder))
  d <- data.table::rbindlist(lapply(f, data.table::fread), fill = TRUE)
  data.table::fwrite(d, cache)
  d
}

# Replicates per (study, method type): expected from the task tables and found on disk.
sim_progress <- function(folder = NULL, root = here::here()) {
  A <- sim_arrays()
  if (!is.null(folder)) A <- A[vapply(A, function(a) a$folder %in% folder, TRUE)]
  do.call(rbind, lapply(A, function(a) {
    tk <- sim_task_table(a$scen, a$designs, a$reps)
    dir <- file.path(root, a$folder, "output", a$study, a$type)
    data.frame(array = a$name, folder = a$folder, tasks = nrow(tk), expected = sum(tk$rep_to - tk$rep_from + 1L),
               found = length(list.files(dir, pattern = "\\.csv$")))
  }))
}

# Design-level means: FDP (with Monte Carlo standard error), power, rejections, incongruity.
design_means <- function(d, by = c("target", "K", "t", "method", "design", "axis", "x")) {
  d <- data.table::as.data.table(d)
  d[, .(FDP = mean(FDP), FDP_se = stats::sd(FDP) / sqrt(.N), n = .N, power = mean(power, na.rm = TRUE),
        R = mean(R), incon_dec = mean(incon_dec, na.rm = TRUE), incon_score = mean(incon_score, na.rm = TRUE)),
    by = by]
}

# Averages over designs: mean FDP, largest design-level FDP, mean power.
overall_means <- function(D, by = c("target", "K", "t", "method")) {
  D <- data.table::as.data.table(D)
  D[, .(FDP = mean(FDP), FDPmax = max(FDP), power = mean(power, na.rm = TRUE), ndes = .N,
        incon_dec = mean(incon_dec, na.rm = TRUE), incon_score = mean(incon_score, na.rm = TRUE)), by = by]
}

# Designs 1, 5, 7 of each target (first and fifth effect-size designs; tau1 = 0.02 for
# pleiotropy and replication, 0.04 for mediation), as in Table 2.
three_designs <- function(d) {
  d <- data.table::as.data.table(d)
  d[design %in% c("effect_01", "effect_05") | (scen != "med" & design == "tau1_002") | (scen == "med" & design == "tau1_004")]
}

# Number formatting for LaTeX tables.
f3 <- function(x, digits = 3) ifelse(is.na(x), "--", formatC(x, format = "f", digits = digits))
fcomma <- function(x) ifelse(is.na(x), "--", gsub(",", "{,}", formatC(x, format = "d", big.mark = ",")))
tex_escape <- function(x) gsub("_", "\\\\_", x)
write_tex <- function(lines, file) { writeLines(lines, file); message("wrote ", file) }

# Table 2 / Supplementary Tables S8-S9: mean FDP (or power) under alternative estimators and
# data-generating processes, averaged over designs 1, 5, 7 (rho = 0.8 block: fifth effect-size
# design of the correlation study). Returns one row per (target, K, t), one column per setting.
SOURCE_ROWS <- list(
  c("pen", "Correlated csmGmm (penalized)"),
  c(NA, "\\textit{Estimation}"),
  c("ML", "Correlated csmGmm, maximum likelihood"),
  c("csm", "csmGmm ($\\bSigma=I$), maximum likelihood"),
  c("csm_cor", "csmGmm, $\\bSigma$ fixed at sample correlation"),
  c("csm_pen", "csmGmm ($\\bSigma=I$), penalized"),
  c("split", "Correlated csmGmm, independent test set"),
  c(NA, "\\textit{Data-generating process}"),
  c("sym", "Zero logistic intercept"),
  c("quant", "Quantitative outcomes"),
  c("exch", "Exchangeable correlation 0.4"),
  c("ld", "Linkage disequilibrium (variant level)"),
  c(NA, "\\textit{Strong correlation ($\\rho=0.8$, fifth effect-size design)}"),
  c("rho08", "Correlated csmGmm"),
  c("rho08_oracle", "Correlated csmGmm, $\\bSigma$ fixed at null correlation"),
  c("rho08_csm", "csmGmm ($\\bSigma=I$), maximum likelihood"))
sources_means <- function(value = "FDP") {
  avg <- function(d, m) {
    if (is.null(d)) return(NULL)
    d <- data.table::as.data.table(d)[method == m]
    if (!nrow(d)) return(NULL)
    g <- d[, .(v = mean(get(value), na.rm = TRUE)), by = .(target, K, t, design)]
    g[, .(v = mean(v)), by = .(target, K, t)]
  }
  opt <- function(folder, study) tryCatch(read_sim(folder, study), error = function(e) { message(conditionMessage(e)); NULL })
  M <- read_sim("Fig2", "main"); S <- read_sim("Fig2", "sym"); C <- opt("Tab2", "main")
  R <- opt("Fig3", "rho"); O <- opt("Tab2", "rho")
  r08 <- function(d) if (is.null(d)) NULL else d[design == "rho_080"]
  td <- function(d) if (is.null(d)) NULL else three_designs(d)
  parts <- list(pen = avg(td(M), "penalized"), ML = avg(td(M), "ML"),
                csm = avg(td(C), "csmGmm"), csm_cor = avg(td(C), "csmGmm_cor"), csm_pen = avg(td(C), "csmGmm_pen"),
                split = avg(td(opt("Tab2", "split")), "penalized"), sym = avg(td(S), "penalized"),
                quant = avg(td(opt("Tab2", "quant")), "penalized"), exch = avg(td(opt("Tab2", "exch")), "penalized"),
                ld = avg(td(opt("Tab2", "ld")), "penalized"),
                rho08 = avg(r08(R), "penalized"), rho08_oracle = avg(r08(O), "penalized_nullSigma"),
                rho08_csm = avg(r08(R), "csmGmm"))
  out <- unique(M[, .(target, K, t)])
  for (nm in names(parts)) {
    p <- parts[[nm]]
    if (is.null(p)) { out[[nm]] <- NA_real_; next }
    data.table::setnames(p, "v", nm); out <- merge(out, p, by = c("target", "K", "t"), all.x = TRUE)
  }
  data.table::setorder(out, target, K, t)
  out
}

# LaTeX version of Table 2 (t = 1) and Supplementary Table S8 (t = 2).
sources_tex <- function(T, t0, label, caption) {
  L <- c("\\begin{table}[!htbp]", "\\centering", sprintf("\\caption{%s}\\label{%s}", caption, label), "\\small",
         "\\begin{tabular}{l cc cc cc}", "\\toprule",
         " & \\multicolumn{2}{c}{Pleiotropy} & \\multicolumn{2}{c}{Replication} & \\multicolumn{2}{c}{Mediation}\\\\",
         "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}\\cmidrule(lr){6-7}",
         " & $K=4$ & $K=6$ & $K=4$ & $K=6$ & $K=4$ & $K=6$\\\\", "\\midrule")
  for (r in SOURCE_ROWS) {
    if (is.na(r[1])) { L <- c(L, sprintf("\\multicolumn{7}{l}{%s}\\\\", r[2])); next }
    cells <- unlist(lapply(c("pleiotropy", "replication", "mediation"), function(tg) sapply(c(4, 6), function(K0) {
      v <- T[target == tg & K == K0 & t == t0][[r[1]]]; f3(if (length(v)) v else NA) })))
    L <- c(L, paste0(if (r[1] == "pen") "" else "\\quad ", r[2], " & ", paste(cells, collapse = " & "), "\\\\"))
  }
  c(L, "\\bottomrule", "\\end{tabular}", "\\end{table}")
}
