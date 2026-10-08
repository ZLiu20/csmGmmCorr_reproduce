# Mediation application with G = 3 mediators (the model of main-text Section 2.1, K = 2G = 6):
#   exposure = a genetic variant (SNP), M1, M2, M3 = the liver expression of three genes near it,
#   outcome Y = LDL cholesterol. One feature per variant, coordinates in mediation order
#   (alpha1, beta1, alpha2, beta2, alpha3, beta3):
#     alpha_g : SNP -> M_g, GTEx v8 liver cis-eQTL z of the variant for gene g (slope / se, ALT allele)
#     beta_g  : M_g -> Y,   S-PrediXcan z of gene g (GTEx v8 MASHR liver model, GLGC EUR LDL)
#   With G = 3, t = 1 and t = 2 are both mediation targets (for t >= G the mediation and
#   pleiotropy nulls coincide; main-text Section 2.1).
#   The genes are chosen from positions only: among the genes with a liver model that were tested
#   with the variant, ordered by the distance of their transcription start site to the variant,
#   rule "nearest" takes the first three; rule "noshare" (the primary data set) takes the nearest
#   gene and then, in order, each gene whose model shares no variant with the models already
#   chosen (with a shared model variant two S-PrediXcan statistics are partly the same statistic).
#   The variants are chosen from positions only: common SNVs (GTEx MAF >= 0.05) with at least three
#   such genes; in each 10-kb window the variant nearest to the middle of the window. No statistic
#   is used to choose features.
#   LDL: GLGC 2021 EUR female and male results combined by inverse-variance weighting; variants
#   are matched to GTEx on the build-37 position and both alleles (GTEx lookup table).
#   S-PrediXcan: z_g = sum_l w_l s_l z_l / s_g, s_l^2 = Gamma_ll, s_g^2 = w' Gamma w over the model
#   variants with an LDL statistic (Gamma = the covariance file of the MASHR model).
# Usage (from the repository root, after bash Data/download_mediation_ldl.sh):
#   Rscript Data/prepare_mediation_ldl.R [<folder of the GTEx and MASHR files>] [<folder that contains the GLGC LDL files>]
#   defaults: Data/raw/mediation_ldl and Data/raw/glgc (searched recursively)
# Output: Data/mediation/ldl_twas_liver.tsv.gz, Data/mediation/mediation_ldl_liver_{noshare,nearest}.tsv.gz,
#   intermediate files in Data/prepared/mediation_ldl/ (resumable: finished steps are skipped).
here::i_am("Data/prepare_mediation_ldl.R")
suppressMessages(library(data.table))
say <- function(...) { cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n"); flush.console() }
a <- commandArgs(trailingOnly = TRUE)
src  <- if (length(a) >= 1) a[1] else here::here("Data", "raw", "mediation_ldl")
groot <- if (length(a) >= 2) a[2] else here::here("Data", "raw", "glgc")
root <- Sys.getenv("MEDLDL_OUT", here::here("Data"))    # tests/test_prepare_mediation_ldl.R writes elsewhere
tmp <- file.path(root, "prepared", "mediation_ldl"); dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
odir <- file.path(root, "mediation"); dir.create(odir, recursive = TRUE, showWarnings = FALSE)
f_pairs <- file.path(src, "Liver.allpairs.txt.gz"); f_look <- file.path(src, "GTEx_v8_lookup_table.txt.gz")
f_db <- file.path(src, "mashr_Liver.db"); f_cov <- file.path(src, "mashr_Liver.txt.gz")
for (f in c(f_pairs, f_look, f_db, f_cov)) if (!file.exists(f)) stop("missing input: ", f)
find1 <- function(pat) { x <- list.files(groot, pattern = pat, recursive = TRUE, full.names = TRUE); if (!length(x)) stop("not found under ", groot, ": ", pat); x[1] }
f_ldl <- c(female = find1("^LDL_INV_EUR_HRC_1KGP3_others_FEMALE\\.meta\\.singlevar\\.results\\.gz$"),
           male   = find1("^LDL_INV_EUR_HRC_1KGP3_others_MALE\\.meta\\.singlevar\\.results\\.gz$"))
zcat <- function(f) paste("gzip -cd", shQuote(f))
strip <- function(g) sub("\\..*$", "", g)
say("inputs: %s ; LDL: %s", src, paste(f_ldl, collapse = ", "))

# ---- 1. liver MASHR model --------------------------------------------------------------------
read_table <- function(db, tab) {
  if (requireNamespace("RSQLite", quietly = TRUE)) {
    con <- DBI::dbConnect(RSQLite::SQLite(), db); on.exit(DBI::dbDisconnect(con))
    as.data.table(DBI::dbReadTable(con, tab))
  } else fread(cmd = sprintf("sqlite3 -header -csv %s 'select * from %s'", shQuote(db), tab))
}
W <- read_table(f_db, "weights"); E <- read_table(f_db, "extra")
say("MASHR weights: %d rows, columns %s", nrow(W), paste(names(W), collapse = ","))
say("MASHR extra: %d genes, columns %s", nrow(E), paste(names(E), collapse = ","))
vcol <- if ("varID" %in% names(W)) "varID" else "rsid"
W <- W[, .(gene = as.character(gene), var = as.character(get(vcol)), eff = toupper(eff_allele), w = as.numeric(weight))]
if (!all(grepl("^chr[0-9XY]+_[0-9]+_[ACGT]+_[ACGT]+_b38$", head(W$var, 1000))))
  stop("unexpected model variant IDs, e.g. ", paste(head(W$var, 3), collapse = " "))
E <- E[, .(gene = as.character(gene), gene_name = as.character(genename))]
W[, c("ref38", "alt38") := tstrsplit(var, "_", keep = 3:4)]

# ---- 2. LDL statistics of the model variants (GLGC EUR, female and male combined) -----------
f_gw <- file.path(tmp, "ldl_model_variants.tsv.gz")
if (!file.exists(f_gw)) {
  vars <- unique(W$var)
  look <- fread(cmd = zcat(f_look), select = c("variant_id", "variant_id_b37"))
  look <- look[variant_id %chin% vars & !is.na(variant_id_b37) & nzchar(variant_id_b37)]
  look[, c("chr37", "pos37", "ref37", "alt37") := tstrsplit(variant_id_b37, "_", keep = 1:4)]
  look[, key := paste(chr37, pos37, sep = ":")]
  say("model variants: %d, with a build-37 position: %d", length(vars), nrow(look))
  rd <- function(f) {
    x <- fread(cmd = zcat(f), select = c("CHROM", "POS_b37", "REF", "ALT", "EFFECT_SIZE", "SE"), colClasses = list(character = c("CHROM", "REF", "ALT")))
    x[, key := paste(sub("^chr", "", CHROM), as.integer(POS_b37), sep = ":")]
    x <- x[key %chin% look$key & is.finite(as.numeric(EFFECT_SIZE)) & as.numeric(SE) > 0]
    x[, .(key, REF = toupper(REF), ALT = toupper(ALT), b = as.numeric(EFFECT_SIZE), s = as.numeric(SE))]
  }
  L <- merge(rd(f_ldl["female"]), rd(f_ldl["male"]), by = c("key", "REF", "ALT"), all = TRUE, suffixes = c("_f", "_m"))
  L[, `:=`(wf = fifelse(is.na(b_f), 0, 1 / s_f^2), wm = fifelse(is.na(b_m), 0, 1 / s_m^2))]
  L[, z := (fifelse(is.na(b_f), 0, b_f * wf) + fifelse(is.na(b_m), 0, b_m * wm)) / sqrt(wf + wm)]
  m <- merge(look[, .(var = variant_id, key, ref37, alt37)], L[, .(key, REF, ALT, z)], by = "key", allow.cartesian = TRUE)
  m <- m[(REF == ref37 & ALT == alt37) | (REF == alt37 & ALT == ref37)]
  m[, z_alt := fifelse(REF == ref37, z, -z)]
  m <- m[!(var %chin% var[duplicated(var)])]       # drop variants matched twice
  fwrite(m[, .(var, z_alt)], f_gw, sep = "\t")
  say("model variants with an LDL statistic: %d", nrow(m))
}
gw <- fread(f_gw)

# ---- 3. S-PrediXcan z of every liver model gene -------------------------------------------
f_tw <- file.path(odir, "ldl_twas_liver.tsv.gz")
if (!file.exists(f_tw)) {
  cv <- fread(cmd = zcat(f_cov)); setnames(cv, toupper(names(cv)))
  if (!all(c("GENE", "RSID1", "RSID2", "VALUE") %in% names(cv))) stop("unexpected covariance columns: ", paste(names(cv), collapse = ","))
  cv <- split(cv[, .(GENE, RSID1, RSID2, VALUE)], by = "GENE", keep.by = FALSE)
  M <- merge(W, gw, by = "var")
  M[, z := fifelse(eff == alt38, z_alt, fifelse(eff == ref38, -z_alt, NA_real_))]
  M <- M[is.finite(z)]
  nmod <- W[, .N, by = gene]
  res <- rbindlist(lapply(split(M, by = "gene"), function(x) {
    g <- x$gene[1]; C <- cv[[g]]; if (is.null(C)) return(NULL)
    v <- x$var; S <- matrix(0, length(v), length(v), dimnames = list(v, v))
    C <- C[RSID1 %chin% v & RSID2 %chin% v]
    S[cbind(C$RSID1, C$RSID2)] <- C$VALUE; S[cbind(C$RSID2, C$RSID1)] <- C$VALUE
    sl <- sqrt(pmax(diag(S), 0)); sg2 <- drop(crossprod(x$w, S %*% x$w))
    if (!(sg2 > 0)) return(NULL)
    data.table(gene = g, n_used = length(v), z_twas = sum(x$w * sl * x$z) / sqrt(sg2))
  }))
  if (nrow(res) < 0.5 * uniqueN(M$gene)) { print(head(rbindlist(cv[1:3], idcol = "GENE"))); print(head(M)); stop("S-PrediXcan failed for most genes: check the variant IDs of the covariance file") }
  res <- merge(res, nmod[, .(gene, n_model = N)], by = "gene")
  res <- merge(res, E, by = "gene", all.x = TRUE)
  res[, gene_id := strip(gene)]
  res <- res[!(gene_id %chin% gene_id[duplicated(gene_id)])]
  fwrite(res[, .(gene_id, gene_name, n_model, n_used, z_twas)], f_tw, sep = "\t")
  say("S-PrediXcan: %d genes; IQR-based SD %.2f; |z| > 5: %d", nrow(res),
      diff(quantile(res$z_twas, c(.25, .75))) / 1.349, sum(abs(res$z_twas) > 5))
  print(head(res[order(-abs(z_twas)), .(gene_name, n_model, n_used, z_twas = round(z_twas, 1))], 25))
}
tw <- fread(f_tw)
modset <- split(W[, .(g = strip(gene), var)], by = "g", keep.by = FALSE)
modset <- lapply(modset, `[[`, "var")

# ---- 4. liver eQTL statistics: one pass over the all-pairs file, split by chromosome --------
sdir <- file.path(tmp, "liver_pairs")
if (!file.exists(file.path(sdir, "DONE"))) {
  dir.create(sdir, showWarnings = FALSE)
  hdr <- names(fread(cmd = paste(zcat(f_pairs), "| head -2")))
  need <- c("gene_id", "variant_id", "tss_distance", "maf", "slope", "slope_se")
  if (!all(need %in% hdr)) stop("unexpected all-pairs columns: ", paste(hdr, collapse = ","))
  ix <- match(need, hdr)
  gfile <- file.path(tmp, "model_genes.txt"); writeLines(tw$gene_id, gfile)
  awk <- sprintf(paste0("awk -F'\\t' -v OFS='\\t' -v G=%s -v D=%s -v cg=%d -v cv=%d -v ct=%d -v cm=%d -v cs=%d -v ce=%d '",
    "BEGIN{while((getline l < G)>0) K[l]=1} NR==1{next} ($cm+0)<0.05{next} ($ce+0)<=0{next} ",
    "{g=$cg; sub(/\\..*$/,\"\",g); if(!(g in K)) next; n=split($cv,a,\"_\"); if(n<5||length(a[3])!=1||length(a[4])!=1) next; ",
    "c=a[1]; sub(/^chr/,\"\",c); if(c !~ /^[0-9]+$/) next; f=D \"/chr\" c \".tsv.gz\"; ",
    "p=\"gzip > \" f; P[p]=1; print g,$cv,a[2],$ct,$cs/$ce | p} END{for(p in P) close(p)}'"), shQuote(gfile), shQuote(sdir), ix[1], ix[2], ix[3], ix[4], ix[5], ix[6])
  say("splitting the all-pairs file (one pass)")
  st <- system(paste(zcat(f_pairs), "|", awk))
  if (st != 0) stop("awk failed")
  writeLines("DONE", file.path(sdir, "DONE"))
}

# ---- 5. features: one variant per 10-kb window, G = 3 genes per variant -------------------
G <- 3L
zt <- setNames(tw$z_twas, tw$gene_id); gn <- setNames(tw$gene_name, tw$gene_id)
nshared <- function(a, b) length(intersect(modset[[a]], modset[[b]]))
pick <- function(genes, rule) {               # genes ordered by distance of the TSS to the variant
  if (rule == "nearest") return(genes[seq_len(G)])
  sel <- genes[1L]
  for (g in genes[-1L]) {
    if (all(vapply(sel, function(s) nshared(s, g) == 0L, TRUE))) sel <- c(sel, g)
    if (length(sel) == G) return(sel)
  }
  NULL
}
out <- list(noshare = list(), nearest = list())
for (k in 1:22) {
  f <- file.path(sdir, sprintf("chr%d.tsv.gz", k)); if (!file.exists(f)) next
  d <- fread(cmd = zcat(f), header = FALSE, col.names = c("gene", "variant", "pos", "tss", "z"))
  d <- d[is.finite(z)]
  d[, ng := .N, by = variant]
  v <- unique(d[ng >= G, .(variant, pos)])
  v[, win := pos %/% 10000L][, dm := abs(pos - (win * 10000L + 5000L))]
  setorder(v, win, dm, pos, variant); v <- v[!duplicated(win)]
  e <- d[variant %chin% v$variant][, ad := abs(tss)]
  setorder(e, variant, ad, gene)
  e <- split(e, by = "variant")
  for (r in names(out)) {
    out[[r]][[k]] <- rbindlist(lapply(e, function(x) {
      g <- pick(x$gene, r); if (is.null(g)) return(NULL)
      m <- match(g, x$gene)
      sh <- sum(combn(g, 2L, function(p) nshared(p[1L], p[2L])))
      row <- list(variant = x$variant[1L], chr = k, pos_b38 = x$pos[1L], shared = sh)
      for (j in seq_len(G)) {
        row[[paste0("gene", j)]] <- g[j]; row[[paste0("gene", j, "_name")]] <- gn[[g[j]]]; row[[paste0("tss", j)]] <- x$tss[m[j]]
        row[[paste0("alpha", j)]] <- x$z[m[j]]; row[[paste0("beta", j)]] <- zt[[g[j]]]
      }
      row
    }))
  }
  say("chr%d: %d variant-gene pairs, %d windows", k, nrow(d), nrow(v))
}
zc <- as.vector(rbind(paste0("alpha", seq_len(G)), paste0("beta", seq_len(G))))
for (r in names(out)) {
  D <- rbindlist(out[[r]]); setorder(D, chr, pos_b38)
  setcolorder(D, c("variant", "chr", "pos_b38", as.vector(rbind(paste0("gene", 1:G), paste0("gene", 1:G, "_name"), paste0("tss", 1:G))), zc, "shared"))
  fwrite(D, file.path(odir, sprintf("mediation_ldl_liver_%s.tsv.gz", r)), sep = "\t")
  Z <- as.matrix(D[, zc, with = FALSE])
  say("rule %s: %d features, %d distinct gene triples, %d with shared model variants", r, nrow(D),
      uniqueN(D[, paste0("gene", 1:G), with = FALSE]), sum(D$shared > 0))
  print(round(cor(Z), 3))
  print(round(apply(Z, 2, function(v) diff(quantile(v, c(.25, .75))) / 1.349), 2))
  print(round(colMeans(abs(Z) > 4), 4))
}
say("DONE")
