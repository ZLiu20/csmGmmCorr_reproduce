# Job arrays of the simulation studies.
#
# Each array runs one method type for one target scenario (pl = pleiotropy, rp = replication,
# med = mediation) in one study. One array task runs CHUNK = 10 consecutive replicates of one
# (design, K, t) cell; the task index (LSB_JOBINDEX) selects the row of the task table.
# For pleiotropy and replication the data do not depend on t (t = 0 in the table), and one
# fit is evaluated at t = 1 and t = 2; mediation data are generated separately for t = 1, 2.
#
# Studies (data-generating processes, SupportingCode/sim_data.R):
#   main  latent AR(1) correlation 0.2, logistic intercept -1, 12 designs per target
#   sym   as main with intercept 0 (sign-symmetric statistics)
#   rho   fifth effect-size design with rho = 0, 0.2, ..., 0.8
#   quant quantitative outcomes;  exch exchangeable correlation 0.4;
#   ld    genotypes in linkage disequilibrium;  split  independent test data set
# Method types (SupportingCode/run_sim_task.R):
#   ours_<sc>_pen   correlated csmGmm, penalized (lambda_mu = 30, lambda_w = 10)
#   ours_<sc>_ml    correlated csmGmm, maximum likelihood
#   ours_<sc>_pen_A / _C            numbers of magnitude components 1/1 and 2/2
#   ours_<sc>_pen_L<l>A<a>          penalty constants lambda_mu = l, lambda_w = a
#   ours_med_pen_oracle             Sigma fixed at the correlation of the null statistics
#   comp_<sc>_csm   csmGmm of Sun et al. (Sigma = I), plug-in Sigma, penalized with Sigma = I
#   comp_pl         QCH and adaFilter;  comp_rp  repfdr and directional partial conjunction
#   comp_med        DEI, MediationTest and the JS-PC combination
CHUNK <- 10L

# designs: rows of design_table(study, scen); reps: replicate range.
sim_task_table <- function(scen, designs, reps = c(1L, 200L), chunk = CHUNK) {
  from <- seq(reps[1], reps[2], by = chunk)
  ch <- data.frame(rep_from = from, rep_to = pmin(reps[2], from + chunk - 1L))
  ts <- if (scen == "med") 1:2 else 0L
  g <- expand.grid(K = c(4L, 6L), design = designs, t = ts, chunk = seq_len(nrow(ch)))
  x <- cbind(g[c("design", "K", "t")], ch[g$chunk, ])
  x <- x[order(-x$K, x$design, x$t, x$rep_from), ]
  rownames(x) <- NULL
  x
}

# Designs 1, 5, 7: first and fifth effect-size designs and the second signal-proportion design.
D3 <- c(1L, 5L, 7L)

# All arrays: name -> folder, study, method type, scenario, designs, replicates, resources.
# Resources: cores, wall time (hh:mm) and memory (GB) per task of 10 replicates.
.arr <- function(folder, study, type, designs, reps, cores = 4L, wall = "12:00", mem = 8L) {
  scen <- regmatches(type, regexpr("(pl|rp|med)", type))
  list(folder = folder, study = study, type = type, scen = scen, designs = designs, reps = reps,
       cores = cores, wall = wall, mem = mem)
}
sim_arrays <- function() {
  A <- list()
  add <- function(name, ...) A[[name]] <<- .arr(...)
  for (sc in c("pl", "rp", "med")) {
    # Figure 2, Table 1, Supplementary Figures 1-12 and Tables 6-7, 16-18
    add(sprintf("main_ours_%s_pen", sc), "Fig2", "main", sprintf("ours_%s_pen", sc), 1:12, c(1L, 500L))
    add(sprintf("main_ours_%s_ml", sc),  "Fig2", "main", sprintf("ours_%s_ml", sc), 1:12, c(1L, 200L))
    add(sprintf("sym_ours_%s_pen", sc),  "Fig2", "sym",  sprintf("ours_%s_pen", sc), 1:12, c(1L, 200L))
    add(sprintf("sym_ours_%s_ml", sc),   "Fig2", "sym",  sprintf("ours_%s_ml", sc), 1:12, c(1L, 200L))
    # Figure 3 and Supplementary Tables S4, S28
    add(sprintf("rho_ours_%s_pen", sc),  "Fig3", "rho",  sprintf("ours_%s_pen", sc), 1:5, c(1L, 200L))
    add(sprintf("rho_ours_%s_ml", sc),   "Fig3", "rho",  sprintf("ours_%s_ml", sc), 1:5, c(1L, 200L))
    add(sprintf("rho_comp_%s_csm", sc),  "Fig3", "rho",  sprintf("comp_%s_csm", sc), 1:5, c(1L, 200L))
    # Table 2 and Supplementary Tables S3, S8-S10
    add(sprintf("main_comp_%s_csm", sc), "Tab2", "main", sprintf("comp_%s_csm", sc), 1:9, c(1L, 200L))
    for (st in c("quant", "exch", "ld", "split"))
      add(sprintf("%s_ours_%s_pen", st, sc), "Tab2", st, sprintf("ours_%s_pen", sc), D3, c(1L, 100L))
    # Supplementary Table S5 (penalty constants) and Table S7 (magnitude components)
    for (pg in c("L0A10", "L10A10", "L100A10", "L30A0", "L30A3", "L30A30"))
      add(sprintf("main_ours_%s_pen_%s", sc, pg), "SuppTab3", "main", sprintf("ours_%s_pen_%s", sc, pg), D3, c(1L, 50L))
    for (s0 in c("A", "C"))
      add(sprintf("main_ours_%s_pen_%s", sc, s0), "SuppTab9", "main", sprintf("ours_%s_pen_%s", sc, s0), D3, c(1L, 50L),
          wall = if (s0 == "C") "16:00" else "12:00")
  }
  add("main_comp_pl",  "Fig2", "main", "comp_pl", 1:12, c(1L, 200L), wall = "4:00", mem = 4L)
  add("main_comp_med", "Fig2", "main", "comp_med", 1:12, c(1L, 200L), cores = 1L, wall = "16:00", mem = 16L)
  add("main_comp_rp",  "SuppTab12", "main", "comp_rp", 1:12, c(1L, 200L), cores = 1L, wall = "12:00", mem = 16L)
  add("rho_ours_med_pen_oracle", "Tab2", "rho", "ours_med_pen_oracle", 4:5, c(1L, 50L))
  for (nm in names(A)) A[[nm]]$name <- nm
  A
}

sim_array <- function(name) {
  A <- sim_arrays()
  if (!name %in% names(A)) stop("unknown array: ", name)
  a <- A[[name]]
  a$tasks <- sim_task_table(a$scen, a$designs, a$reps)
  a
}
