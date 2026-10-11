CHUNK <- 10L

sim_task_table <- function(scen, designs, reps = c(1L, 200L), chunk = CHUNK) {
  starts <- seq.int(as.integer(reps[1L]), as.integer(reps[2L]), by = as.integer(chunk))
  ranges <- data.frame(
    rep_from = starts,
    rep_to = pmin(as.integer(reps[2L]), starts + as.integer(chunk) - 1L)
  )
  thresholds <- if (identical(scen, "med")) 1:2 else 0L
  grid <- expand.grid(
    K = c(4L, 6L),
    design = as.integer(designs),
    t = thresholds,
    block = seq_len(nrow(ranges)),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  out <- cbind(grid[c("design", "K", "t")], ranges[grid$block, ])
  out <- out[order(out$K, out$design, out$t, out$rep_from), , drop = FALSE]
  rownames(out) <- NULL
  out
}

.key_designs <- c(1L, 5L, 7L)

.array_spec <- function(folder, study, type, designs, reps,
                        cores = 4L, wall = "12:00", memory = 8L) {
  scenario <- regmatches(type, regexpr("(pl|rp|med)", type))
  list(
    folder = folder,
    study = study,
    type = type,
    scen = scenario,
    designs = as.integer(designs),
    reps = as.integer(reps),
    cores = as.integer(cores),
    wall = wall,
    mem = as.integer(memory)
  )
}

sim_arrays <- function() {
  arrays <- list()
  register <- function(name, ...) arrays[[name]] <<- .array_spec(...)

  for (scenario in c("pl", "rp", "med")) {
    register(sprintf("main_ours_%s_pen", scenario), "Fig2", "main",
             sprintf("ours_%s_pen", scenario), 1:12, c(1L, 500L))
    register(sprintf("main_ours_%s_ml", scenario), "Fig2", "main",
             sprintf("ours_%s_ml", scenario), 1:12, c(1L, 200L))
    register(sprintf("sym_ours_%s_pen", scenario), "Fig2", "sym",
             sprintf("ours_%s_pen", scenario), 1:12, c(1L, 200L))
    register(sprintf("sym_ours_%s_ml", scenario), "Fig2", "sym",
             sprintf("ours_%s_ml", scenario), 1:12, c(1L, 200L))
    register(sprintf("rho_ours_%s_pen", scenario), "Fig3", "rho",
             sprintf("ours_%s_pen", scenario), 1:5, c(1L, 200L))
    register(sprintf("rho_ours_%s_ml", scenario), "Fig3", "rho",
             sprintf("ours_%s_ml", scenario), 1:5, c(1L, 200L))
    register(sprintf("rho_comp_%s_csm", scenario), "Fig3", "rho",
             sprintf("comp_%s_csm", scenario), 1:5, c(1L, 200L))
    register(sprintf("main_comp_%s_csm", scenario), "Tab2", "main",
             sprintf("comp_%s_csm", scenario), 1:9, c(1L, 200L))

    for (study in c("quant", "exch", "ld", "split")) {
      register(sprintf("%s_ours_%s_pen", study, scenario), "Tab2", study,
               sprintf("ours_%s_pen", scenario), .key_designs, c(1L, 100L))
    }

    for (penalty in c("L0A10", "L10A10", "L100A10", "L30A0", "L30A3", "L30A30")) {
      register(sprintf("main_ours_%s_pen_%s", scenario, penalty), "SuppTab3", "main",
               sprintf("ours_%s_pen_%s", scenario, penalty), .key_designs, c(1L, 50L))
    }

    for (shape in c("A", "C")) {
      register(sprintf("main_ours_%s_pen_%s", scenario, shape), "SuppTab9", "main",
               sprintf("ours_%s_pen_%s", scenario, shape), .key_designs, c(1L, 50L),
               wall = if (shape == "C") "16:00" else "12:00")
    }
  }

  register("main_comp_pl", "Fig2", "main", "comp_pl", 1:12, c(1L, 500L),
           wall = "4:00", memory = 4L)
  register("main_comp_med", "Fig2", "main", "comp_med", 1:12, c(1L, 200L),
           cores = 1L, wall = "16:00", memory = 16L)
  register("main_comp_rp", "SuppTab12", "main", "comp_rp", 1:12, c(1L, 200L),
           cores = 1L, wall = "12:00", memory = 16L)
  register("rho_ours_med_pen_oracle", "Tab2", "rho", "ours_med_pen_oracle", 4:5,
           c(1L, 50L))

  for (name in names(arrays)) arrays[[name]]$name <- name
  arrays
}

sim_array <- function(name) {
  arrays <- sim_arrays()
  if (!name %in% names(arrays)) stop("Unknown simulation array: ", name)
  ans <- arrays[[name]]
  ans$tasks <- sim_task_table(ans$scen, ans$designs, ans$reps)
  ans
}
