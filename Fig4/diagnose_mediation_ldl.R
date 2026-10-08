# Diagnostics of the LDL mediation fits (Fig4/output/mediation_ldl/*), t = 1:
#   rejected variants; those with a strong pair (|alpha_g| >= 4 and |beta_g| >= 4 for some g);
#   those with all |alpha_g| < 2; loci (greedy 500-kb clumping by lfdr) and loci that contain at
#   least one rejected variant with a strong pair. The last count asks whether discoveries without
#   a strong pair are linkage-disequilibrium neighbours of a strong pathway (locus-level reading)
#   or stand alone. Reference: strong-pair variants and their loci in the data.
#   Rscript Fig4/diagnose_mediation_ldl.R
here::i_am("Fig4/diagnose_mediation_ldl.R")
source(here::here("SupportingCode", "load_all.R"))
suppressMessages(library(data.table))
strong <- function(d) with(d, (abs(alpha1) >= 4 & abs(beta1) >= 4) | (abs(alpha2) >= 4 & abs(beta2) >= 4) |
                              (abs(alpha3) >= 4 & abs(beta3) >= 4))
# locus of every selected variant: the nearest lead on its chromosome (all are within 500 kb)
loci_of <- function(chr, pos, score, win = 5e5) {
  lead <- clump_leads(chr, pos, score, rep(TRUE, length(pos)), win)
  id <- integer(length(pos))
  for (c0 in unique(chr)) {
    i <- which(chr == c0); l <- lead[chr[lead] == c0]
    id[i] <- l[vapply(pos[i], function(p) which.min(abs(pos[l] - p)), 1L)]
  }
  id
}
od <- here::here("Fig4", "output", "mediation_ldl")
res <- rbindlist(lapply(list.dirs(od, recursive = FALSE), function(v) {
  f <- file.path(v, "discoveries.tsv.gz"); if (!file.exists(f)) return(NULL)
  d <- fread(f)[t == 1 & rej == TRUE]
  if (!nrow(d)) return(data.table(fit = basename(v), rejected = 0L))
  s <- strong(d); ma <- with(d, pmax(abs(alpha1), abs(alpha2), abs(alpha3)))
  loc <- loci_of(as.character(d$chr), d$pos_b38, d$lfdr)
  data.table(fit = basename(v), rejected = nrow(d), strong_pair = sum(s), all_alpha_lt2 = sum(ma < 2),
             loci = uniqueN(loc), loci_with_strong = uniqueN(loc[s]),
             variants_in_strong_loci = sum(loc %in% loc[s]))
}), fill = TRUE)
print(res, nrows = 100)
for (set in c("noshare", "noshare_dedup")) {
  x <- freadgz(here::here("Data", "mediation", "mediation_ldl_liver_noshare.tsv.gz"))
  if (set == "noshare_dedup") x <- medldl_dedup(x)
  s <- strong(x)
  cat(sprintf("data %s: %d variants, %d with a strong pair, in %d loci\n", set, nrow(x), sum(s),
              length(clump_leads(as.character(x$chr[s]), x$pos_b38[s], -abs(x$beta1[s]), rep(TRUE, sum(s))))))
}
