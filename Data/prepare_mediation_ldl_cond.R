# Conditional exposure-mediator statistics for the LDL mediation application.
# The marginal liver eQTL z of a variant j for gene g picks up r * z_c from the lead eQTL variant c
# of g when j and c are in linkage disequilibrium (r = their correlation), so variants without an
# effect of their own have nonzero marginal alpha. Here alpha is replaced by the z-statistic of j
# conditional on the lead variant of g (one-step conditioning, as in COJO):
#   alpha_cond = (z_j - r z_c) / sqrt(1 - r^2),
# with r from the 1000 Genomes EUR reference. Exceptions: j = c keeps its marginal z; when
# r^2 >= R2_PROXY, j is a proxy of the lead and keeps its marginal z (the two cannot be told apart).
# Features whose r with any lead is unavailable (variant or lead not in the reference) are dropped.
# The lead variant of a gene is its common SNV with the largest |z| in the liver eQTL file.
# beta (S-PrediXcan) is unchanged.
# Input: Data/mediation/mediation_ldl_liver_noshare.tsv.gz and Data/prepared/mediation_ldl/liver_pairs/
# (Data/prepare_mediation_ldl.R), the GTEx v8 lookup table (build-37 positions), and a PLINK 1.9
# binary fileset of the 1000 Genomes EUR samples on build 37.
#   Rscript Data/prepare_mediation_ldl_cond.R <1000G EUR bfile prefix> [plink binary]
# Output: Data/mediation/mediation_ldl_liver_noshare_cond.tsv.gz (alpha1-3 conditional; the
# marginal values in alpha1_marg-alpha3_marg; how each alpha was obtained in alpha1_how-alpha3_how).
here::i_am("Data/prepare_mediation_ldl_cond.R")
suppressMessages(library(data.table))
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n"); flush.console() }
a <- commandArgs(trailingOnly = TRUE)
bfile <- a[1]; plink <- if (length(a) >= 2) a[2] else "plink"
if (!file.exists(paste0(bfile, ".bim"))) stop("PLINK fileset not found: ", bfile)
R2_PROXY <- 0.8
zcat <- function(f) paste("gzip -cd", shQuote(f))
root <- Sys.getenv("MEDLDL_OUT", here::here("Data"))
tmp <- file.path(root, "prepared", "mediation_ldl"); sdir <- file.path(tmp, "liver_pairs")
cdir <- file.path(tmp, "cond"); dir.create(cdir, recursive = TRUE, showWarnings = FALSE)
src <- Sys.getenv("MEDLDL_SRC", here::here("Data", "raw", "mediation_ldl"))
D <- fread(cmd = zcat(file.path(root, "mediation", "mediation_ldl_liver_noshare.tsv.gz")))
genes <- unique(c(D$gene1, D$gene2, D$gene3))
say("features: %d, genes: %d", nrow(D), length(genes))

# ---- 1. lead eQTL variant of each gene --------------------------------------------------------
leads <- rbindlist(lapply(1:22, function(k) {
  f <- file.path(sdir, sprintf("chr%d.tsv.gz", k)); if (!file.exists(f)) return(NULL)
  d <- fread(cmd = zcat(f), header = FALSE, col.names = c("gene", "variant", "pos", "tss", "z"))
  d <- d[gene %chin% genes & is.finite(z)]
  d[d[, .I[which.max(abs(z))], by = gene]$V1, .(gene, lead = variant, lead_z = z)]
}))
say("lead variants: %d genes", nrow(leads))

# ---- 2. build-37 positions and alleles (GTEx lookup table) ------------------------------------
vars <- unique(c(D$variant, leads$lead))
look <- fread(cmd = zcat(file.path(src, "GTEx_v8_lookup_table.txt.gz")), select = c("variant_id", "variant_id_b37"))
look[, v := sub("_b38$", "", variant_id)]
look <- look[v %chin% vars & !is.na(variant_id_b37) & nzchar(variant_id_b37)]
look[, c("c37", "p37", "r37", "a37") := tstrsplit(variant_id_b37, "_", keep = 1:4)]
look <- look[!(v %chin% v[duplicated(v)])]
say("variants: %d, with a build-37 position: %d", length(vars), nrow(look))

# ---- 3. subset of the reference, variant IDs = GTEx IDs, allele signs --------------------------
bim <- fread(paste0(bfile, ".bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "a1", "a2"),
             colClasses = list(character = c(1, 2, 5, 6)))
bim[, chr := sub("^chr", "", chr)][, key := paste(chr, pos, sep = ":")]
look[, key := paste(c37, p37, sep = ":")]
m <- merge(look[, .(v, key, r37, a37)], bim[, .(key, chr, pos, a1, a2)], by = "key", allow.cartesian = TRUE)
m <- m[(a1 == a37 & a2 == r37) | (a1 == r37 & a2 == a37)]
m[, s := fifelse(a1 == a37, 1, -1)]              # +1 when the reference A1 allele is the GTEx ALT allele
m <- m[!(v %chin% v[duplicated(v)]) & !(key %chin% key[duplicated(key)])]
say("variants matched to the reference (position and alleles): %d", nrow(m))
rng <- file.path(cdir, "range.txt"); fwrite(m[, .(chr, pos, pos, v)], rng, sep = " ", col.names = FALSE)
sub_b <- file.path(cdir, "ref_subset")
st <- system2(plink, c("--bfile", bfile, "--extract", "range", rng, "--make-bed", "--out", sub_b, "--allow-no-sex"))
if (st != 0) stop("plink --extract failed")
sb <- fread(paste0(sub_b, ".bim"), header = FALSE, col.names = c("chr", "id", "cm", "pos", "a1", "a2"),
            colClasses = list(character = c(1, 2, 5, 6)))
sb[, chr := sub("^chr", "", chr)][, key := paste(chr, pos, sep = ":")]
sb[, id := m$v[match(paste(key, a1, a2), paste(m$key, m$a1, m$a2))]]
sb[is.na(id), id := sprintf("unused_%d", .I)]
fwrite(sb[, .(chr, id, cm, pos, a1, a2)], paste0(sub_b, ".bim"), sep = "\t", col.names = FALSE)
lf <- file.path(cdir, "leads.txt"); writeLines(intersect(unique(leads$lead), sb$id), lf)
ldo <- file.path(cdir, "lead_ld")
st <- system2(plink, c("--bfile", sub_b, "--r", "--ld-snp-list", lf, "--ld-window-kb", "1500", "--ld-window", "999999",
                       "--ld-window-r2", "0", "--out", ldo, "--allow-no-sex"))
if (st != 0) stop("plink --r failed")
ld <- fread(paste0(ldo, ".ld"), select = c("SNP_A", "SNP_B", "R"))
sg <- setNames(m$s, m$v)
ld[, R := R * sg[SNP_A] * sg[SNP_B]]           # correlation of the ALT-allele counts
ld <- ld[is.finite(R)]
ld <- unique(rbind(ld[, .(lead = SNP_A, variant = SNP_B, r = R)], ld[, .(lead = SNP_B, variant = SNP_A, r = R)]),
             by = c("lead", "variant"))
setkey(ld, lead, variant)
say("lead-variant correlations: %d pairs", nrow(ld))

# ---- 4. conditional alpha -----------------------------------------------------------------------
lz <- setNames(leads$lead_z, leads$gene); lv <- setNames(leads$lead, leads$gene)
cond_alpha <- function(variant, gene, z) {
  cl <- unname(lv[gene]); zc <- unname(lz[gene])
  r <- ld[.(cl, variant), r]
  how <- fifelse(!is.na(cl) & variant == cl, "lead", fifelse(is.na(r), "no_ld", fifelse(r^2 >= R2_PROXY, "proxy", "conditional")))
  val <- fifelse(how %chin% c("lead", "proxy"), z, (z - r * zc) / sqrt(1 - r^2))
  val[how == "no_ld"] <- NA_real_
  list(val, how, r)
}
for (g in 1:3) {
  res <- cond_alpha(D$variant, D[[paste0("gene", g)]], D[[paste0("alpha", g)]])
  D[, (paste0("alpha", g, "_marg")) := get(paste0("alpha", g))]
  D[, (paste0("alpha", g)) := res[[1]]][, (paste0("alpha", g, "_how")) := res[[2]]][, (paste0("r_lead", g)) := res[[3]]]
}
# check of the allele alignment: among conditioned variants with |r| > 0.3 the marginal z must
# follow r * z_lead (a negative correlation would mean a sign error)
chk <- rbindlist(lapply(1:3, function(g) D[get(paste0("alpha", g, "_how")) == "conditional" & abs(get(paste0("r_lead", g))) > 0.3,
  .(zm = get(paste0("alpha", g, "_marg")), pred = get(paste0("r_lead", g)) * lz[get(paste0("gene", g))])]))
cc <- if (nrow(chk) > 10) cor(chk$zm, chk$pred) else NA
say("alignment check: cor(marginal z, r * z_lead) = %.3f over %d values", cc, nrow(chk))
if (is.na(cc) || cc < 0.3) stop("allele alignment check failed")
print(rbindlist(lapply(1:3, function(g) D[, .N, by = .(alpha = g, how = get(paste0("alpha", g, "_how")))])))
D <- D[stats::complete.cases(D[, .(alpha1, alpha2, alpha3)])]
fwrite(D, file.path(root, "mediation", "mediation_ldl_liver_noshare_cond.tsv.gz"), sep = "\t")
Zm <- as.matrix(D[, .(alpha1_marg, alpha2_marg, alpha3_marg)]); Zc <- as.matrix(D[, .(alpha1, alpha2, alpha3)])
iq <- function(x) diff(stats::quantile(x, c(.25, .75))) / 1.349
say("features kept: %d", nrow(D))
print(data.table(alpha = 1:3, iqr_sd_marginal = round(apply(Zm, 2, iq), 2), iqr_sd_conditional = round(apply(Zc, 2, iq), 2),
                 abs_gt4_marginal = round(colMeans(abs(Zm) > 4), 4), abs_gt4_conditional = round(colMeans(abs(Zc) > 4), 4)))
say("DONE")
