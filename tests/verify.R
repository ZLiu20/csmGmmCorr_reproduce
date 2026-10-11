# Compare the regenerated tables and figure data with the values of the paper in reference/.
#   Rscript tests/verify.R            (run by the final job of run_all.sh, after make_all.sh)
# reference/ mirrors the repository: reference/Tab1/Tab1.tex is compared with Tab1/Tab1.tex, etc.
#   .tex, .txt   all numbers in the file (thousands separators removed) and the text around them
#   .csv, .tsv   same columns and rows; numeric columns compared value by value, others exactly;
#                columns with run times (minutes, seconds, elapsed, time, date, *_min) are not
#                compared; columns of the reference that the new run does not write are listed
# Run times printed in tables (VOLATILE below) are removed before the comparison.
# Status per file: IDENTICAL (same content), EQUAL (numbers agree within TOL, default 1e-6
# relative), CLOSE (every number within NOISE_REL = 1% relative or NOISE_ABS = 0.005 absolute:
# the size of the run-to-run variation of multithreaded fits, see README), DIFFERENT (with the
# largest discrepancy), MISSING (not regenerated).
# Also compares logs/versions.txt with reference/versions.txt. Report: logs/verify_report.txt
# and logs/verify_report.tsv. The exit status is 1 if a file is DIFFERENT or MISSING.
here::i_am("tests/verify.R")
suppressMessages(library(data.table))
tol <- as.numeric(Sys.getenv("TOL", "1e-6"))
noise_rel <- as.numeric(Sys.getenv("NOISE_REL", "0.01")); noise_abs <- as.numeric(Sys.getenv("NOISE_ABS", "0.005"))
root <- here::here(); ref_root <- file.path(root, "reference")
if (!dir.exists(ref_root)) stop("reference/ not found")
files <- list.files(ref_root, recursive = TRUE, all.files = FALSE)
files <- setdiff(files, c("versions.txt", "README.md", "raw_md5.txt", "files.txt"))
dir.create(file.path(root, "logs"), showWarnings = FALSE)
# file -> regular expression of run times to remove before comparing (last columns "Iter." and "Min." of the
# table of all real-data fits)
VOLATILE <- list("SuppTab19/SuppTab20.tex" = "& [0-9{},]+ & [0-9]+\\.[0-9](?=\\\\\\\\$)")
TIMECOL <- "minute|second|elapsed|time|date|_min$|^iter"
num_re <- "[-+]?(?:[0-9]+\\.?[0-9]*|\\.[0-9]+)(?:[eE][-+]?[0-9]+)?"
clean <- function(x) gsub("\\{,\\}|(?<=[0-9]),(?=[0-9]{3})", "", x, perl = TRUE)
numbers <- function(x) as.numeric(regmatches(x, gregexpr(num_re, x, perl = TRUE))[[1]])
reldiff <- function(a, b) { d <- abs(a - b); s <- pmax(abs(a), abs(b)); ifelse(d == 0, 0, d / pmax(s, 1e-12)) }
cmp_num <- function(a, b) {
  ok <- is.na(a) == is.na(b); keep <- !is.na(a) & !is.na(b); a <- a[keep]; b <- b[keep]
  if (!all(ok)) return(list(rel = Inf, abs = Inf, n = sum(!ok), far = sum(!ok)))
  r <- reldiff(a, b); far <- sum(r > noise_rel & abs(a - b) > noise_abs)
  list(rel = if (length(r)) max(r) else 0, abs = if (length(r)) max(abs(a - b)) else 0, n = sum(r > tol), far = far)
}
compare_text <- function(fr, fn, volatile = NULL) {
  rd <- function(f) { x <- readLines(f, warn = FALSE); if (!is.null(volatile)) x <- gsub(volatile, "& ", x, perl = TRUE); clean(paste(x, collapse = "\n")) }
  a <- rd(fr); b <- rd(fn)
  if (identical(a, b)) return(list(status = "IDENTICAL", detail = ""))
  ta <- gsub(num_re, "#", a, perl = TRUE); tb <- gsub(num_re, "#", b, perl = TRUE)
  na <- numbers(a); nb <- numbers(b)
  if (!identical(ta, tb) || length(na) != length(nb))
    return(list(status = "DIFFERENT", detail = sprintf("text or layout differs (%d vs %d numbers)", length(na), length(nb))))
  r <- cmp_num(na, nb)
  if (r$n == 0) list(status = "EQUAL", detail = sprintf("max rel. diff %.2g", r$rel))
  else if (r$far == 0) list(status = "CLOSE", detail = sprintf("%d of %d numbers differ, all within noise; max rel. %.2g", r$n, length(na), r$rel))
  else { i <- which.max(reldiff(na, nb)); list(status = "DIFFERENT", detail = sprintf("%d of %d numbers differ; largest: %s vs %s", r$n, length(na), format(na[i]), format(nb[i]))) }
}
compare_table <- function(fr, fn) {
  if (identical(unname(tools::md5sum(fr)), unname(tools::md5sum(fn)))) return(list(status = "IDENTICAL", detail = ""))
  a <- fread(fr); b <- fread(fn)
  if (nrow(a) != nrow(b)) return(list(status = "DIFFERENT", detail = sprintf("%d vs %d rows", nrow(a), nrow(b))))
  cols <- grep(TIMECOL, intersect(names(a), names(b)), ignore.case = TRUE, value = TRUE, invert = TRUE)
  absent <- grep(TIMECOL, setdiff(names(a), names(b)), ignore.case = TRUE, value = TRUE, invert = TRUE)
  worst <- 0; bad <- character(0); far <- 0
  for (v in cols) {
    x <- a[[v]]; y <- b[[v]]
    if (is.numeric(x) && is.numeric(y)) { r <- cmp_num(as.numeric(x), as.numeric(y)); worst <- max(worst, r$rel); far <- far + r$far; if (r$n) bad <- c(bad, sprintf("%s (%d rows, max rel. %.2g)", v, r$n, r$rel)) }
    else if (!identical(as.character(x), as.character(y))) { far <- far + 1; bad <- c(bad, sprintf("%s (%d rows)", v, sum(as.character(x) != as.character(y), na.rm = TRUE))) }
  }
  note <- if (length(absent)) sprintf("; reference columns not written: %s", paste(absent, collapse = ", ")) else ""
  if (!length(bad)) list(status = "EQUAL", detail = sprintf("max rel. diff %.2g%s", worst, note))
  else list(status = if (far == 0) "CLOSE" else "DIFFERENT", detail = paste0(paste(head(bad, 4), collapse = "; "), note))
}
res <- rbindlist(lapply(files, function(f) {
  fr <- file.path(ref_root, f); fn <- file.path(root, f)
  o <- if (!file.exists(fn)) list(status = "MISSING", detail = "")
       else if (grepl("\\.(csv|tsv)(\\.gz)?$", f)) tryCatch(compare_table(fr, fn), error = function(e) list(status = "DIFFERENT", detail = conditionMessage(e)))
       else compare_text(fr, fn, VOLATILE[[f]])
  data.table(file = f, status = o$status, detail = o$detail)
}))
res[, status := factor(status, c("MISSING", "DIFFERENT", "CLOSE", "EQUAL", "IDENTICAL"))]
setorder(res, status, file)
fwrite(res, file.path(root, "logs", "verify_report.tsv"), sep = "\t")
out <- c(sprintf("Comparison with reference/ (%s), tolerance %g", format(Sys.time(), "%F %T"), tol), "",
         capture.output(print(res[, .N, by = status], row.names = FALSE)), "")
vr <- file.path(ref_root, "versions.txt"); vn <- file.path(root, "logs", "versions.txt")
if (file.exists(vr) && file.exists(vn)) {
  a <- readLines(vr); b <- readLines(vn); d <- setdiff(union(a, b), intersect(a, b))
  out <- c(out, if (length(d)) c("Versions differ from the reference:", paste("  reference:", setdiff(a, b)), paste("  this run: ", setdiff(b, a)))
                else "R and package versions: as in the reference", "")
}
mf <- file.path(root, "logs", "make", "status.txt")
if (file.exists(mf)) out <- c(out, "Scripts (make_all.sh):", readLines(mf), "")
out <- c(out, capture.output(print(res, row.names = FALSE, nrows = 1000)))
writeLines(out, file.path(root, "logs", "verify_report.txt")); writeLines(out)
quit(status = as.integer(any(res$status %in% c("MISSING", "DIFFERENT"))))
