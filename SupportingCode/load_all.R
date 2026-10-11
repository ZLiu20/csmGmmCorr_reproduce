.supporting_files <- c(
  "sim_data.R", "incongruity.R", "sim_tasks.R", "run_sim_task.R",
  "summarize.R", "realdata.R", "glgc_prepare.R", "plot_style.R", "floor_study.R"
)
for (.file in .supporting_files) source(here::here("SupportingCode", .file))

.comparator_files <- list.files(here::here("SupportingCode", "comparators"),
                                pattern = "[.]R$", full.names = TRUE)
for (.file in .comparator_files) source(.file)
rm(.file, .supporting_files, .comparator_files)

load_kernels <- function(root = here::here()) {
  extension <- .Platform$dynlib.ext
  library_path <- file.path(root, "SupportingCode", "src", paste0("reprokernels", extension))
  if (!file.exists(library_path)) stop("Compiled reproduction kernels are missing: ", library_path)
  if (!is.loaded("rd_majorant")) dyn.load(library_path)
  invisible(library_path)
}
