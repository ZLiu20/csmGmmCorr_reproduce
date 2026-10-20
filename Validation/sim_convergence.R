# Convergence of the csmGmm fits in all simulation studies. A replicate whose two starts both fail
# to converge makes no discoveries and enters the reported FDP and power as such (Supplementary
# Section S4.1); this script counts how often that happened.
# Writes Validation/sim_convergence_by_design.csv (one row per study, target, K, t, design and
# method) and prints the overall failure rates and the largest design-level failure rate.
#   Rscript Validation/sim_convergence.R      (after the simulation arrays; uses the cached tables)
here::i_am("Validation/sim_convergence.R")
suppressMessages(library(data.table))
source(here::here("SupportingCode", "load_all.R"))
A <- sim_arrays()
fs <- unique(rbindlist(lapply(A, function(a) data.table(folder = a$folder, study = a$study,
                                                          fit = grepl("^ours|_csm$", a$type)))))
fs <- unique(fs[fit == TRUE, .(folder, study)])
out <- list()
for (i in seq_len(nrow(fs))) {
  d <- tryCatch(read_sim(fs$folder[i], fs$study[i]), error = function(e) NULL)
  if (is.null(d) || !"converged" %in% names(d)) next
  d <- d[!is.na(converged) & !grepl("_adjusted$", method)]
  # one row per fitted replicate: mediation data are generated for each t, the other targets once
  d[, t_data := if ("scen" %in% names(d)) fifelse(scen == "med", t, 0L) else t]
  d <- unique(d, by = c("scen", "design", "K", "rep", "t_data", "method"))
  d[, `:=`(folder = fs$folder[i], study_ = fs$study[i])]
  out[[i]] <- d[, .(n = .N, failed = sum(!as.logical(converged))),
                by = .(folder, study = study_, scen, K, t = t_data, design, method)]
}
res <- rbindlist(out)
fwrite(res, here::here("Validation", "sim_convergence_by_design.csv"))
tot <- res[, .(replicates = sum(n), failed = sum(failed), rate = sum(failed) / sum(n),
               max_design_rate = max(failed / n)), by = .(folder, study, scen, method)][order(-rate)]
options(width = 160)
cat(sprintf("fitted replicates: %d; failed: %d (%.4f%%)\n", sum(res$n), sum(res$failed), 100 * sum(res$failed) / sum(res$n)))
print(tot, nrows = 200)
cat("\ndesigns with any failure:\n"); print(res[failed > 0][order(-failed / n)], nrows = 60)
