# Parse every R file of the repository (syntax check; no code is run).
files <- list.files(".", pattern = "\\.R$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("SupportingCode/packages/", files)]
bad <- 0
for (f in files) { e <- tryCatch({ parse(f); NULL }, error = function(e) conditionMessage(e)); if (!is.null(e)) { bad <- bad + 1; cat("PARSE ERROR", f, "\n", e, "\n") } }
cat(sprintf("%d R files parsed, %d with errors\n", length(files), bad))
if (bad) quit(status = 1)
