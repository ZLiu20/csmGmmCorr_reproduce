# Preparation of the GLGC summary statistics (Supplementary Section S5.2).
#
# For each analysis (manifest Data/glgc_<analysis>_manifest.tsv, one row per coordinate):
#   * per source: build-37 autosomal positions, REF/ALT upper case (ALT is the effect allele),
#     z = beta / se (or, if the standard error is unusable, the signed normal quantile of the
#     p-value), minor allele frequency >= 0.01, N >= half of the largest N of the source;
#     a repeated (chr, pos, REF, ALT) key keeps the first record with the largest N;
#   * keys present in every source (exact chromosome, position, REF and ALT; no strand
#     flipping and no sign reversal);
#   * statistics capped at +-12;
#   * one fitting variant per occupied 10-kb window: the variant whose key has the smallest
#     BLAKE2b-64 hash (deterministic; independent of the statistics).
# Output in <out_dir>: chr<k>.tsv.gz (variant_id, rsid, chr, pos_b37, ref, alt, one column per
# source), fit_sample.tsv.gz, chromosome_manifest.tsv, source_qc.tsv, prepare_summary.tsv.
# The files are identical (after decompression) to those of the streaming Python
# implementation used for the manuscript, which wrote CRLF line ends.
GLGC_MISSING <- c("", "NA", "NaN", "nan", ".")

# Acklam's inverse normal approximation (used by the original implementation; |error| < 1e-9).
acklam_qnorm <- function(p) {
  a <- c(-3.969683028665376e01, 2.209460984245205e02, -2.759285104469687e02, 1.383577518672690e02,
         -3.066479806614716e01, 2.506628277459239e00)
  b <- c(-5.447609879822406e01, 1.615858368580409e02, -1.556989798598866e02, 6.680131188771972e01,
         -1.328068155288572e01)
  cc <- c(-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e00, -2.549732539343734e00,
          4.374664141464968e00, 2.938163982698783e00)
  d <- c(7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e00, 3.754408661907416e00)
  out <- rep(NaN, length(p)); low <- 0.02425; high <- 1 - low
  tail_fun <- function(q) ((((((cc[1] * q + cc[2]) * q + cc[3]) * q + cc[4]) * q + cc[5]) * q + cc[6]) /
                           ((((d[1] * q + d[2]) * q + d[3]) * q + d[4]) * q + 1))
  i <- which(p > 0 & p < low);   if (length(i)) out[i] <- tail_fun(sqrt(-2 * log(p[i])))
  i <- which(p > high & p < 1);  if (length(i)) out[i] <- -tail_fun(sqrt(-2 * log(1 - p[i])))
  i <- which(p >= low & p <= high)
  if (length(i)) {
    q <- p[i] - 0.5; r <- q * q
    out[i] <- (((((a[1] * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * r + a[6]) * q /
              (((((b[1] * r + b[2]) * r + b[3]) * r + b[4]) * r + b[5]) * r + 1)
  }
  out[p == 0] <- -Inf; out[p == 1] <- Inf
  out
}

glgc_chrom_rank <- function(x) {
  v <- toupper(sub("^chr", "", trimws(x), ignore.case = TRUE))
  r <- rep(10000L, length(v))
  num <- grepl("^[+-]?[0-9]+$", v)
  r[num] <- suppressWarnings(as.integer(v[num]))
  r[v == "X"] <- 23L; r[v == "Y"] <- 24L; r[v %in% c("M", "MT")] <- 25L
  r[is.na(r)] <- 10000L
  r
}
glgc_chrom_label <- function(rank) ifelse(rank == 23L, "X", ifelse(rank == 24L, "Y", ifelse(rank == 25L, "MT", as.character(rank))))

.num <- function(x) if (is.numeric(x)) as.double(x) else suppressWarnings(as.numeric(ifelse(x %in% GLGC_MISSING, NA, x)))

# One source: QC as described above; returns the retained keys with z and rsID, and counts.
glgc_read_source <- function(path, beta_col, se_col, p_col, af_col, n_col, min_maf, min_n_fraction,
                             max_chrom_rank = 22L) {
  need <- c("rsID", "CHROM", "POS_b37", "REF", "ALT", beta_col, se_col, af_col, n_col)
  hdr <- names(data.table::fread(cmd = paste("gzip -cd", shQuote(path), "| head -1"), sep = "\t"))
  miss <- setdiff(need, hdr); if (length(miss)) stop(path, " missing columns: ", paste(miss, collapse = ", "))
  sel <- unique(c(need, if (p_col %in% hdr) p_col))
  x <- data.table::fread(cmd = paste("gzip -cd", shQuote(path)), sep = "\t", select = sel,
                         colClasses = list(character = c("rsID", "CHROM", "REF", "ALT")),
                         na.strings = GLGC_MISSING, showProgress = FALSE)
  cnt <- list(raw_rows = nrow(x))
  nval <- .num(x[[n_col]])
  max_n <- suppressWarnings(max(nval[is.finite(nval)]))
  if (!is.finite(max_n) || max_n <= 0) stop("no valid sample size in ", path)
  min_n <- max_n * min_n_fraction
  rank <- glgc_chrom_rank(x$CHROM)
  pos <- suppressWarnings(trunc(.num(x$POS_b37)))
  bad <- is.na(pos) | rank < 1L | rank > 25L | pos < 1
  cnt$invalid_coordinate <- sum(bad)
  excl <- !bad & rank > max_chrom_rank
  cnt$non_autosomal_or_excluded <- sum(excl)
  keep <- !bad & !excl
  ref <- toupper(trimws(ifelse(is.na(x$REF), "", x$REF))); alt <- toupper(trimws(ifelse(is.na(x$ALT), "", x$ALT)))
  bad_allele <- keep & (ref == "" | alt == "" | ref == alt)
  cnt$invalid_allele <- sum(bad_allele); keep <- keep & !bad_allele
  beta <- .num(x[[beta_col]]); se <- .num(x[[se_col]])
  pv <- if (p_col %in% names(x)) .num(x[[p_col]]) else rep(NA_real_, nrow(x))
  z <- ifelse(is.finite(beta) & is.finite(se) & se > 0, beta / se, NA_real_)
  usep <- !(is.finite(beta) & is.finite(se) & se > 0) & is.finite(pv) & pv > 0 & pv <= 1 & is.finite(beta)
  if (any(usep & keep)) {
    i <- which(usep & keep); pp <- pmin(pmax(pv[i], 1e-300), 1)
    q <- acklam_qnorm(1 - pp / 2); z[i] <- ifelse(beta[i] < 0, -q, q)
  }
  nz <- keep & !is.finite(z); cnt$nonfinite_z <- sum(nz); keep <- keep & !nz
  af <- .num(x[[af_col]])
  maf_fail <- keep & !(is.finite(af) & af >= min_maf & af <= 1 - min_maf)
  cnt$maf_fail <- sum(maf_fail); keep <- keep & !maf_fail
  n_fail <- keep & !(is.finite(nval) & nval >= min_n)
  cnt$n_fail <- sum(n_fail); keep <- keep & !n_fail
  d <- data.table::data.table(rank = rank[keep], pos = as.integer(pos[keep]), ref = ref[keep], alt = alt[keep],
                              z = z[keep], n = nval[keep], rsid = x$rsID[keep], idx = which(keep))
  rm(x); gc(verbose = FALSE)
  # Repeated keys: the streaming implementation kept the first record and replaced it only by
  # a record with strictly larger N, i.e. it kept the first record with the largest N.
  data.table::setorderv(d, c("rank", "pos", "ref", "alt", "idx"))
  d[, cm := cummax(n), by = .(rank, pos, ref, alt)]
  d[, prev := data.table::shift(cm, fill = -Inf), by = .(rank, pos, ref, alt)]
  dup <- d[, .N, by = .(rank, pos, ref, alt)][N > 1L]
  cnt$duplicate_allele_at_position <- sum(dup$N - 1L)
  cnt$qc_rows <- sum(d$n > d$prev)            # accepted records, including replacements
  d <- d[d[, .I[which.max(n)], by = .(rank, pos, ref, alt)]$V1]   # which.max: first maximum
  d[, c("cm", "prev", "n", "idx") := NULL]
  list(data = d, counts = cnt, max_n = max_n, min_n = min_n)
}

prepare_glgc <- function(manifest, raw_root, out_dir, min_maf = 0.01, min_n_fraction = 0.5,
                         fit_window_bp = 10000L, z_cap = 12, overwrite = FALSE) {
  suppressMessages(library(data.table))
  load_kernels()
  man <- fread(manifest, sep = "\t", colClasses = "character")
  if (nrow(man) < 2L || anyDuplicated(man$label)) stop("need at least two sources with unique labels")
  if (file.exists(file.path(out_dir, "PREPARATION_COMPLETE")) && !overwrite)
    stop("completed output exists: ", out_dir, " (use overwrite = TRUE)")
  unlink(out_dir, recursive = TRUE); dir.create(out_dir, recursive = TRUE)
  D <- NULL; qc <- list()
  for (i in seq_len(nrow(man))) {
    m <- man[i]; path <- if (startsWith(m$relative_path, "/")) m$relative_path else file.path(raw_root, m$relative_path)
    message(sprintf("%s  %-11s reading %s", format(Sys.time(), "%H:%M:%S"), m$label, basename(path)))
    s <- glgc_read_source(path, m$beta_col, m$se_col, if (nzchar(m$p_col)) m$p_col else "pvalue",
                          if (nzchar(m$af_col)) m$af_col else "POOLED_ALT_AF", if (nzchar(m$n_col)) m$n_col else "N",
                          min_maf, min_n_fraction)
    qc[[i]] <- data.table(label = m$label, path = normalizePath(path), max_n = s$max_n, min_n = s$min_n,
                          min_maf = min_maf, as.data.table(s$counts))
    d <- s$data; setnames(d, c("z", "rsid"), c(m$label, paste0("rsid_", i)))
    D <- if (is.null(D)) d else merge(D, d, by = c("rank", "pos", "ref", "alt"))
    message(sprintf("           %d records pass QC; %d keys common so far", nrow(s$data), nrow(D)))
    rm(s, d); gc(verbose = FALSE)
  }
  setorderv(D, c("rank", "pos", "ref", "alt"))   # C-locale order of the alleles, as in Python
  rs <- rep(".", nrow(D))
  for (i in rev(seq_len(nrow(man)))) {          # first non-missing rsID in manifest order
    v <- D[[paste0("rsid_", i)]]; ok <- !is.na(v) & !(v %in% GLGC_MISSING); rs[ok] <- v[ok]
  }
  lab <- glgc_chrom_label(D$rank)
  out <- data.table(variant_id = sprintf("chr%s:%d:%s:%s", lab, D$pos, D$ref, D$alt), rsid = rs, chr = lab,
                    pos_b37 = as.character(D$pos), ref = D$ref, alt = D$alt)
  for (l in man$label) out[[l]] <- sprintf("%.10g", pmin(pmax(D[[l]], -z_cap), z_cap))
  # fitting sample: smallest hash per window (rank, (pos - 1) %/% fit_window_bp)
  h <- .Call("glgc_blake2b8_hex", out$variant_id)
  win <- paste(D$rank, (D$pos - 1L) %/% as.integer(fit_window_bp))
  o <- order(win, h, method = "radix")
  first <- o[!duplicated(win[o])]
  fit_rows <- sort(first)                        # output order
  wr <- function(x, f) fwrite(x, f, sep = "\t", eol = "\r\n", quote = "auto", compress = if (grepl("gz$", f)) "gzip" else "none")
  wr(out[fit_rows], file.path(out_dir, "fit_sample.tsv.gz"))
  chm <- NULL
  for (r in sort(unique(D$rank))) {
    sel <- which(D$rank == r); f <- sprintf("chr%s.tsv.gz", glgc_chrom_label(r))
    wr(out[sel], file.path(out_dir, f))
    chm <- rbind(chm, data.table(chr = glgc_chrom_label(r), path = f, n_rows = length(sel)))
  }
  wr(chm, file.path(out_dir, "chromosome_manifest.tsv"))
  wr(rbindlist(qc), file.path(out_dir, "source_qc.tsv"))
  wr(data.table(strict_common = nrow(out), fit_sample_n = length(fit_rows), fit_window_bp = fit_window_bp,
                min_maf = min_maf, min_n_fraction = min_n_fraction, z_cap = z_cap, n_sources = nrow(man),
                manifest = normalizePath(manifest), raw_root = normalizePath(raw_root)),
     file.path(out_dir, "prepare_summary.tsv"))
  writeLines("ok", file.path(out_dir, "PREPARATION_COMPLETE"))
  message(sprintf("PREPARATION_COMPLETE common = %d, fit_sample = %d, out = %s", nrow(out), length(fit_rows), out_dir))
  invisible(list(common = nrow(out), fit = length(fit_rows)))
}
