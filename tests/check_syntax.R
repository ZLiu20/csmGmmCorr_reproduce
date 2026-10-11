root <- normalizePath(if (length(commandArgs(trailingOnly = TRUE))) commandArgs(trailingOnly = TRUE)[1L] else ".")
rfiles <- list.files(root, pattern = "[.]R$", recursive = TRUE, full.names = TRUE)
for (file in rfiles) tryCatch(parse(file), error = function(e) stop(file, ": ", conditionMessage(e)))
arrays <- new.env(parent = baseenv())
sys.source(file.path(root, "SupportingCode", "sim_tasks.R"), arrays)
a <- arrays$sim_arrays()
stopifnot(length(a) == 64L, length(unique(vapply(a, `[[`, "", "name"))) == length(a))
wrappers <- list.files(root, pattern = "_sim[.]R$", recursive = TRUE, full.names = TRUE)
stopifnot(length(wrappers) >= 64L)

for (name in names(a)) {
  definition <- a[[name]]
  wrapper_rel <- file.path(definition$folder, paste0(name, "_sim.R"))
  lsf_rel <- file.path(definition$folder, paste0("run_", name, ".lsf"))
  wrapper_path <- file.path(root, wrapper_rel)
  lsf_path <- file.path(root, lsf_rel)
  stopifnot(file.exists(wrapper_path), file.exists(lsf_path))

  wrapper_text <- readLines(wrapper_path, warn = FALSE)
  stopifnot(
    any(grepl(sprintf('run_sim_task("%s"', name), wrapper_text, fixed = TRUE)),
    any(grepl(wrapper_rel, wrapper_text, fixed = TRUE))
  )

  task_count <- nrow(arrays$sim_array(name)$tasks)
  lsf_text <- readLines(lsf_path, warn = FALSE)
  stopifnot(
    any(grepl(sprintf("[1-%d]", task_count), lsf_text, fixed = TRUE)),
    any(grepl(sprintf("#BSUB -W %s", definition$wall), lsf_text, fixed = TRUE)),
    any(grepl(sprintf("#BSUB -n %d", definition$cores), lsf_text, fixed = TRUE)),
    any(grepl(sprintf("#BSUB -M %d", definition$mem), lsf_text, fixed = TRUE)),
    any(grepl(wrapper_rel, lsf_text, fixed = TRUE))
  )
}

shell_files <- list.files(root, pattern = "[.](sh|lsf)$", recursive = TRUE,
                          full.names = TRUE)
bash <- Sys.which("bash")
if (nzchar(bash)) {
  for (file in shell_files) {
    status <- system2(bash, c("-n", file), stdout = FALSE, stderr = FALSE)
    if (status != 0L) stop("Shell syntax failed: ", file)
  }
}
source(file.path(root, "SupportingCode", "realdata.R"))
stopifnot(nrow(RD_JOBS) == 13L, nrow(MEDLDL_JOBS) == 23L,
          length(unique(basename(rd_fit_dirs("output")))) == 13L,
          length(unique(basename(medldl_fit_dirs("output")))) == 23L)
run_all <- readLines(file.path(root, "run_all.sh"), warn = FALSE)
stopifnot(any(grepl("p\\$array\\[p\\$found >= p\\$expected\\]", run_all)),
          any(grepl("count_realdata_fits", run_all, fixed = TRUE)))
cat(sprintf("SYNTAX_AND_LAYOUT_OK R=%d arrays=%d wrappers=%d shell=%d\n",
            length(rfiles), length(a), length(wrappers), length(shell_files)))
