GLGC_MISSING <- c("", "NA", "NaN", "nan", ".")

# The statistical fallback uses R's high-accuracy normal quantile.  It is
# applied only when a valid beta and two-sided p-value are available but the
# reported standard error cannot produce a finite z-statistic.
acklam_qnorm <- function(p) stats::qnorm(p)

glgc_chrom_rank <- function(x) {
  value <- toupper(sub("^CHR", "", trimws(as.character(x)), ignore.case = TRUE))
  out <- suppressWarnings(as.integer(value))
  out[value == "X"] <- 23L
  out[value == "Y"] <- 24L
  out[value %in% c("M", "MT")] <- 25L
  out[is.na(out)] <- 10000L
  out
}

glgc_chrom_label <- function(rank) {
  ifelse(rank == 23L, "X",
         ifelse(rank == 24L, "Y",
                ifelse(rank == 25L, "MT", as.character(rank))))
}

.glgc_num <- function(x) {
  if (is.numeric(x)) return(as.double(x))
  suppressWarnings(as.numeric(ifelse(x %in% GLGC_MISSING, NA, x)))
}

.glgc_hash <- function(keys) {
  if (!is.loaded("glgc_blake2b8_hex")) {
    library_path <- here::here(
      "SupportingCode", "src", paste0("reprokernels", .Platform$dynlib.ext)
    )
    if (!file.exists(library_path)) {
      stop("Compiled GLGC hash helper is missing; run Rscript setup.R first.")
    }
    dyn.load(library_path)
  }
  .Call("glgc_blake2b8_hex", as.character(keys))
}

glgc_read_source <- function(path, beta_col, se_col, p_col, af_col, n_col,
                             min_maf, min_n_fraction,
                             max_chrom_rank = 22L) {
  command <- paste("gzip -cd", shQuote(path))
  header <- names(data.table::fread(cmd = command, nrows = 0L, sep = "\t"))
  required <- c("rsID", "CHROM", "POS_b37", "REF", "ALT",
                beta_col, se_col, af_col, n_col)
  missing <- setdiff(required, header)
  if (length(missing)) {
    stop(path, " is missing columns: ", paste(missing, collapse = ", "))
  }
  selected <- unique(c(required, if (nzchar(p_col)) p_col))
  x <- data.table::fread(
    cmd = command, sep = "\t", select = selected,
    colClasses = list(character = c("rsID", "CHROM", "REF", "ALT")),
    na.strings = GLGC_MISSING, showProgress = FALSE
  )

  n_value <- .glgc_num(x[[n_col]])
  finite_n <- n_value[is.finite(n_value)]
  if (!length(finite_n)) stop("No valid sample size in ", path)
  max_n <- max(finite_n)
  min_n <- max_n * min_n_fraction
  rank <- glgc_chrom_rank(x$CHROM)
  position <- suppressWarnings(as.integer(x$POS_b37))
  ref <- toupper(trimws(ifelse(is.na(x$REF), "", x$REF)))
  alt <- toupper(trimws(ifelse(is.na(x$ALT), "", x$ALT)))

  valid_coordinate <- is.finite(position) & position > 0L
  autosomal <- rank >= 1L & rank <= max_chrom_rank
  valid_allele <- nzchar(ref) & nzchar(alt) & ref != alt
  beta <- .glgc_num(x[[beta_col]])
  se <- .glgc_num(x[[se_col]])
  p_value <- if (nzchar(p_col) && p_col %in% names(x)) {
    .glgc_num(x[[p_col]])
  } else {
    rep(NA_real_, nrow(x))
  }
  z <- ifelse(is.finite(beta) & is.finite(se) & se > 0, beta / se, NA_real_)
  use_p <- !is.finite(z) & is.finite(p_value) & p_value > 0 &
    p_value <= 1 & is.finite(beta)
  if (any(use_p)) {
    tail_probability <- pmin(pmax(p_value[use_p], 1e-300), 1)
    quantile <- acklam_qnorm(1 - tail_probability / 2)
    z[use_p] <- ifelse(beta[use_p] < 0, -quantile, quantile)
  }
  af <- .glgc_num(x[[af_col]])
  maf <- pmin(af, 1 - af)
  valid_maf <- is.finite(maf) & maf >= min_maf
  valid_n <- is.finite(n_value) & n_value >= min_n

  base <- valid_coordinate & autosomal & valid_allele
  keep <- base & is.finite(z) & valid_maf & valid_n
  d <- data.table::data.table(
    rank = rank[keep], pos = position[keep], ref = ref[keep], alt = alt[keep],
    z = z[keep], n = n_value[keep], rsid = x$rsID[keep],
    input_order = which(keep)
  )

  data.table::setorder(d, rank, pos, ref, alt, -n, input_order)
  duplicate_rows <- d[, .N, by = .(rank, pos, ref, alt)][N > 1L, sum(N - 1L)]
  d <- d[, .SD[1L], by = .(rank, pos, ref, alt)]
  d[, c("n", "input_order") := NULL]

  counts <- data.table::data.table(
    raw_rows = nrow(x),
    invalid_coordinate = sum(!valid_coordinate),
    non_autosomal_or_excluded = sum(valid_coordinate & !autosomal),
    invalid_allele = sum(valid_coordinate & autosomal & !valid_allele),
    invalid_z = sum(base & !is.finite(z)),
    maf_fail = sum(base & is.finite(z) & !valid_maf),
    n_fail = sum(base & is.finite(z) & valid_maf & !valid_n),
    duplicate_key_rows = ifelse(length(duplicate_rows), duplicate_rows, 0L),
    retained = nrow(d)
  )
  list(data = d, counts = counts, max_n = max_n, min_n = min_n)
}

prepare_glgc <- function(manifest, raw_root, out_dir, min_maf = 0.01,
                         min_n_fraction = 0.5, fit_window_bp = 10000L,
                         z_cap = 12, overwrite = FALSE) {
  man <- data.table::fread(manifest, sep = "\t", colClasses = "character")
  if (nrow(man) < 2L || anyDuplicated(man$label)) {
    stop("Manifest requires at least two uniquely labeled sources.")
  }
  marker <- file.path(out_dir, "PREPARATION_COMPLETE")
  if (file.exists(marker) && !overwrite) {
    stop("Prepared output already exists; use overwrite=TRUE to replace it.")
  }
  unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  common <- NULL
  qc <- vector("list", nrow(man))
  for (i in seq_len(nrow(man))) {
    item <- man[i]
    path <- if (startsWith(item$relative_path, "/")) {
      item$relative_path
    } else {
      file.path(raw_root, item$relative_path)
    }
    message(sprintf("[%s] reading %s", item$label, basename(path)))
    source <- glgc_read_source(
      path, item$beta_col, item$se_col, item$p_col,
      if (nzchar(item$af_col)) item$af_col else "POOLED_ALT_AF",
      if (nzchar(item$n_col)) item$n_col else "N",
      min_maf, min_n_fraction
    )
    qc[[i]] <- cbind(
      data.table::data.table(label = item$label,
                             path = normalizePath(path, mustWork = TRUE),
                             max_n = source$max_n, min_n = source$min_n),
      source$counts
    )
    d <- data.table::copy(source$data)
    data.table::setnames(d, c("z", "rsid"),
                         c(item$label, paste0("rsid_", i)))
    common <- if (is.null(common)) d else {
      merge(common, d, by = c("rank", "pos", "ref", "alt"))
    }
  }

  data.table::setorder(common, rank, pos, ref, alt)
  rsid_columns <- paste0("rsid_", seq_len(nrow(man)))
  rsid <- rep(NA_character_, nrow(common))
  for (column in rev(rsid_columns)) {
    value <- common[[column]]
    usable <- !is.na(value) & nzchar(value) & !value %in% GLGC_MISSING
    rsid[usable] <- value[usable]
  }
  chromosome <- glgc_chrom_label(common$rank)
  out <- data.table::data.table(
    variant_id = sprintf("chr%s:%d:%s:%s", chromosome, common$pos,
                         common$ref, common$alt),
    rsid = rsid, chr = chromosome, pos_b37 = common$pos,
    ref = common$ref, alt = common$alt
  )
  for (label in man$label) {
    out[[label]] <- sprintf("%.10g", pmin(pmax(common[[label]], -z_cap), z_cap))
  }

  window <- paste(common$rank, common$pos %/% as.integer(fit_window_bp))
  hash <- .glgc_hash(out$variant_id)
  ordered <- order(window, hash, method = "radix")
  selected <- sort(ordered[!duplicated(window[ordered])])

  write_table <- function(x, path) {
    data.table::fwrite(x, path, sep = "\t", quote = FALSE, eol = "\r\n")
  }
  write_table(out[selected], file.path(out_dir, "fit_sample.tsv.gz"))
  chromosome_manifest <- lapply(sort(unique(common$rank)), function(rank) {
    rows <- which(common$rank == rank)
    label <- glgc_chrom_label(rank)
    filename <- sprintf("chr%s.tsv.gz", label)
    write_table(out[rows], file.path(out_dir, filename))
    data.table::data.table(chr = label, path = filename, n_rows = length(rows))
  })
  write_table(data.table::rbindlist(chromosome_manifest),
              file.path(out_dir, "chromosome_manifest.tsv"))
  write_table(data.table::rbindlist(qc, fill = TRUE),
              file.path(out_dir, "source_qc.tsv"))
  write_table(
    data.table::data.table(
      strict_common = nrow(out), fit_sample = length(selected),
      fit_window_bp = fit_window_bp, min_maf = min_maf,
      min_n_fraction = min_n_fraction, z_cap = z_cap,
      n_sources = nrow(man), manifest = normalizePath(manifest),
      raw_root = normalizePath(raw_root)
    ),
    file.path(out_dir, "prepare_summary.tsv")
  )
  writeLines("ok", marker)
  message(sprintf("PREPARATION_COMPLETE common=%d fit_sample=%d out=%s",
                  nrow(out), length(selected), out_dir))
  invisible(list(common = nrow(out), fit = length(selected)))
}
