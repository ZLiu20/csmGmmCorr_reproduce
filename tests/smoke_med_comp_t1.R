# Smoke test: mediation comparators at t = 1 with and without the partial-conjunction screen.
# For a few designs of the main study, the same data are analysed by DEI, MediationTest and JS
# with the earlier combination max(partial conjunction, pathway Bonferroni) ("old") and with the
# pathway Bonferroni p-value alone ("new", now the default at t = 1).
# Checks: (1) new p-values <= old p-values for every feature; (2) the new BH rejection set
# contains the old one; (3) the FDP of the new rule. Small N and few replicates: for testing only.
#   SMOKE_N=20000 SMOKE_REPS=2 Rscript tests/smoke_med_comp_t1.R
here::i_am("tests/smoke_med_comp_t1.R")
source(here::here("SupportingCode", "load_all.R"))
N <- as.integer(Sys.getenv("SMOKE_N", "20000")); REPS <- seq_len(as.integer(Sys.getenv("SMOKE_REPS", "2")))
DESIGNS <- as.integer(strsplit(Sys.getenv("SMOKE_DESIGNS", "1,5,9,12"), ",")[[1]])
Q <- 0.1; t <- 1L
FNS <- list(DEI = run_dei_mediation_R1, MediationTest = run_mediation_test_R1, JS = run_mediation_safeguard_R1)
out <- list()
for (K in c(4L, 6L)) for (d in DESIGNS) for (r in REPS) {
  dat <- make_dataset("med", d, K, t, r, N, study = "main")$test
  alt <- truth_alt(dat$H, "mediation", t)
  for (m in names(FNS)) {
    options(csm.med_pc_screen_t1 = TRUE);  o_old <- FNS[[m]](dat$Z, t, Q)
    options(csm.med_pc_screen_t1 = FALSE); o_new <- FNS[[m]](dat$Z, t, Q)
    fdp <- function(rj) if (sum(rj)) mean(!alt[rj]) else 0
    pw <- function(rj) if (sum(alt)) sum(rj & alt) / sum(alt) else NA
    out[[length(out) + 1L]] <- data.frame(K = K, design = d, rep = r, method = m, nAlt = sum(alt),
      R_old = sum(o_old$rejected), R_new = sum(o_new$rejected),
      FDP_old = fdp(o_old$rejected), FDP_new = fdp(o_new$rejected),
      power_old = pw(o_old$rejected), power_new = pw(o_new$rejected),
      p_not_larger = all(o_new$score <= o_old$score + 1e-15),
      superset = all(o_new$rejected[o_old$rejected]),
      tuning_new = o_new$tuning)
    cat(sprintf("K=%d design=%d rep=%d %-13s R %d -> %d  FDP %.3f -> %.3f  power %.3f -> %.3f  checks %s %s\n",
                K, d, r, m, sum(o_old$rejected), sum(o_new$rejected), fdp(o_old$rejected), fdp(o_new$rejected),
                pw(o_old$rejected), pw(o_new$rejected), tail(out, 1)[[1]]$p_not_larger, tail(out, 1)[[1]]$superset))
  }
}
res <- do.call(rbind, out)
dir.create(here::here("tests", "logs"), showWarnings = FALSE)
utils::write.csv(res, here::here("tests", "logs", "smoke_med_comp_t1.csv"), row.names = FALSE)
cat("\nall p-values not larger:", all(res$p_not_larger), "; all rejection sets nested:", all(res$superset), "\n")
print(aggregate(cbind(FDP_old, FDP_new, power_old, power_new) ~ method + K, data = res, FUN = mean), digits = 3)
cat("tuning (new):", unique(res$tuning_new), sep = "\n  ")
