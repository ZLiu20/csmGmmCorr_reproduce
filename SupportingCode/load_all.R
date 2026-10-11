# Sources all supporting code (base R only; packages are loaded by the scripts that need them).
local({
  code <- here::here("SupportingCode")
  for (f in c("sim_data.R", "incongruity.R", "sim_tasks.R", "run_sim_task.R", "summarize.R",
              "realdata.R", "glgc_prepare.R", "plot_style.R", "floor_study.R"))
    source(file.path(code, f))
  for (f in list.files(file.path(code, "comparators"), pattern = "\\.R$", full.names = TRUE)) source(f)
})
# Compiled kernels (monotone adjustment and incongruity for millions of variants; BLAKE2b
# hash for the GLGC fitting sample). Built once by setup.R.
load_kernels <- function() {
  lib <- here::here("SupportingCode", "src", paste0("reprokernels", .Platform$dynlib.ext))
  if (!file.exists(lib)) stop("compiled kernels not found; run  Rscript setup.R  in the repository root")
  if (!is.loaded("rd_majorant")) dyn.load(lib)
  invisible(lib)
}
