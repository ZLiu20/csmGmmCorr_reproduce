# Test of Data/prepare_mediation_ldl.R on a small synthetic data set (runs in seconds).
#   Rscript tests/test_prepare_mediation_ldl.R
# Checks the allele matching of the LDL statistics, the S-PrediXcan formula, the choice of the
# variant per 10-kb window and of the two genes (rules "nearest" and "noshare").
here::i_am("tests/test_prepare_mediation_ldl.R")
suppressMessages(library(data.table))
if (!requireNamespace("RSQLite", quietly = TRUE)) stop("this test needs the R package RSQLite (install.packages('RSQLite'))")
td <- file.path(tempdir(), "medldl_test"); unlink(td, recursive = TRUE)
src <- file.path(td, "src"); gl <- file.path(td, "glgc"); out <- file.path(td, "out")
dir.create(src, recursive = TRUE); dir.create(gl); dir.create(out)
wgz <- function(x, f, sep = "\t") { fwrite(x, sub("\\.gz$", "", f), sep = sep); system(paste("gzip -f", shQuote(sub("\\.gz$", "", f)))) }

# model: G1 = {v1, v2}, G2 = {v2 (effect allele = REF)}, G3 = {v3}, G4 = {v4}, G5 = {v5}
v <- data.table(var = c("chr1_1000_A_G_b38", "chr1_2000_C_T_b38", "chr1_3000_G_A_b38", "chr2_1000_A_C_b38", "chr2_5000_T_G_b38"),
                b37 = c("1_1100_A_G_b37", "1_2100_C_T_b37", "1_3100_G_A_b37", "2_1100_A_C_b37", "2_5100_T_G_b37"))
W <- data.table(gene = c("ENSG01.1", "ENSG01.1", "ENSG02.1", "ENSG03.1", "ENSG04.1", "ENSG05.1"),
                rsid = v$var[c(1, 2, 2, 3, 4, 5)], varID = v$var[c(1, 2, 2, 3, 4, 5)],
                ref_allele = c("A", "C", "T", "G", "A", "T"), eff_allele = c("G", "T", "C", "A", "C", "G"),
                weight = c(0.5, -0.3, 0.4, 0.7, -0.2, 0.6))
E <- data.table(gene = paste0("ENSG0", 1:5, ".1"), genename = paste0("GENE", 1:5), gene_type = "protein_coding",
                n.snps.in.model = c(2, 1, 1, 1, 1), pred.perf.R2 = NA, pred.perf.pval = NA, pred.perf.qval = NA)
con <- DBI::dbConnect(RSQLite::SQLite(), file.path(src, "mashr_Liver.db"))
DBI::dbWriteTable(con, "weights", as.data.frame(W)); DBI::dbWriteTable(con, "extra", as.data.frame(E)); DBI::dbDisconnect(con)
cv <- data.table(GENE = c("ENSG01.1", "ENSG01.1", "ENSG01.1", "ENSG02.1", "ENSG03.1", "ENSG04.1", "ENSG05.1"),
                 RSID1 = v$var[c(1, 1, 2, 2, 3, 4, 5)], RSID2 = v$var[c(1, 2, 2, 2, 3, 4, 5)], VALUE = c(0.30, 0.10, 0.40, 0.40, 0.25, 0.20, 0.35))
wgz(cv, file.path(src, "mashr_Liver.txt.gz"), sep = " ")
p38 <- tstrsplit(v$var, "_")
wgz(data.table(variant_id = v$var, chr = p38[[1]], variant_pos = p38[[2]], ref = p38[[3]], alt = p38[[4]], num_alt_per_site = 1,
               rs_id_dbSNP151_GRCh38p7 = paste0("rs", 1:5), variant_id_b37 = v$b37), file.path(src, "GTEx_v8_lookup_table.txt.gz"))
# GLGC: v2 is reported with the alleles swapped; v5 only in the male file
mk <- function(i, b, s) { p <- tstrsplit(v$b37[i], "_")
  data.table(CHROM = p[[1]], POS_b37 = as.integer(p[[2]]), rsID = paste0("rs", i), REF = p[[3]], ALT = p[[4]], EFFECT_SIZE = b, SE = s,
             pvalue = 0.5, POOLED_ALT_AF = 0.3, N = 1000) }
fe <- mk(1:4, c(0.02, 0.03, -0.01, 0.05), c(0.01, 0.01, 0.02, 0.01)); fe[2, c("REF", "ALT") := .("T", "C")]
ma <- mk(1:5, c(0.04, 0.01, 0.00, 0.03, -0.06), c(0.02, 0.01, 0.01, 0.01, 0.02)); ma[2, c("REF", "ALT") := .("T", "C")]
wgz(fe, file.path(gl, "LDL_INV_EUR_HRC_1KGP3_others_FEMALE.meta.singlevar.results.gz"))
wgz(ma, file.path(gl, "LDL_INV_EUR_HRC_1KGP3_others_MALE.meta.singlevar.results.gz"))
# eQTL all-pairs: chr1 variants every 500 bp from 500 to 29,500; genes G1 (TSS 5,000), G2 (8,000),
# G3 (25,000) and G9 without a model (6,000); chr2 genes G4 (TSS 1,000) and G5 (4,000)
set.seed(1)
ap <- rbindlist(list(
  CJ(gene_id = c("ENSG01.1", "ENSG02.1", "ENSG03.1", "ENSG09.1"), pos = seq(500, 29500, 500))[, chr := "chr1"],
  CJ(gene_id = c("ENSG04.1", "ENSG05.1"), pos = seq(500, 9500, 500))[, chr := "chr2"]))
tss <- c(ENSG01.1 = 5000, ENSG02.1 = 8000, ENSG03.1 = 25000, ENSG09.1 = 6000, ENSG04.1 = 1000, ENSG05.1 = 4000)
ap[, `:=`(variant_id = sprintf("%s_%d_A_G_b38", chr, pos), tss_distance = pos - tss[gene_id], ma_samples = 10, ma_count = 10,
          maf = ifelse(pos == 15000, 0.01, 0.2), pval_nominal = 0.5, slope = round(rnorm(.N), 3), slope_se = 0.5)]
wgz(ap[, .(gene_id, variant_id, tss_distance, ma_samples, ma_count, maf, pval_nominal, slope, slope_se)], file.path(src, "Liver.allpairs.txt.gz"))

st <- system2("Rscript", c(here::here("Data", "prepare_mediation_ldl.R"), shQuote(src), shQuote(gl)), env = paste0("MEDLDL_OUT=", out))
if (st != 0) stop("prepare_mediation_ldl.R failed")

ok <- function(cond, msg) { cat(if (cond) "OK  " else "FAIL", msg, "\n"); if (!cond) assign("nfail", nfail + 1, .GlobalEnv) }
nfail <- 0
# expected S-PrediXcan z
zf <- fe$EFFECT_SIZE / fe$SE; zm <- ma$EFFECT_SIZE / ma$SE
ivw <- function(i) { b <- c(fe$EFFECT_SIZE[i], ma$EFFECT_SIZE[i]); s <- c(fe$SE[i], ma$SE[i]); ok <- !is.na(b); sum(b[ok] / s[ok]^2) / sqrt(sum(1 / s[ok]^2)) }
zalt <- c(ivw(1), -ivw(2), ivw(3), ivw(4), ma$EFFECT_SIZE[5] / ma$SE[5])
G1 <- { w <- c(0.5, -0.3); S <- matrix(c(0.30, 0.10, 0.10, 0.40), 2); sum(w * sqrt(diag(S)) * zalt[1:2]) / sqrt(drop(t(w) %*% S %*% w)) }
G2 <- sign(0.4) * (-zalt[2])                     # one variant, effect allele = REF
G3 <- sign(0.7) * zalt[3]; G4 <- sign(-0.2) * zalt[4]; G5 <- sign(0.6) * zalt[5]
tw <- fread(file.path(out, "mediation", "ldl_twas_liver.tsv.gz"))
ok(nrow(tw) == 5, "S-PrediXcan: 5 genes")
ok(isTRUE(all.equal(tw[order(gene_id)]$z_twas, c(G1, G2, G3, G4, G5), tolerance = 1e-6)), "S-PrediXcan z (allele flips, IVW, covariance)")
ns <- fread(file.path(out, "mediation", "mediation_ldl_liver_noshare.tsv.gz"))
ne <- fread(file.path(out, "mediation", "mediation_ldl_liver_nearest.tsv.gz"))
# chr1 windows: [0,10k) -> 5,000 ; [10k,20k) -> 15,000 has MAF 0.01, so 14,500 ; [20k,30k) -> 25,000; chr2 [0,10k) -> 5,000
ok(identical(ne$pos_b38, c(5000L, 14500L, 25000L, 5000L)), "one variant per 10-kb window, nearest to the middle, MAF filter")
r1 <- ne[chr == 1 & pos_b38 == 5000]
ok(r1$gene1 == "ENSG01" && r1$gene2 == "ENSG02" && r1$shared == 1, "nearest: G1 then G2 (shared model variant); gene without a model skipped")
r2 <- ns[chr == 1 & pos_b38 == 5000]
ok(r2$gene1 == "ENSG01" && r2$gene2 == "ENSG03" && r2$shared == 0, "noshare: G1 then G3")
r3 <- ne[chr == 1 & pos_b38 == 25000]
ok(r3$gene1 == "ENSG03" && r3$gene2 == "ENSG02", "nearest at 25,000: G3 then G2")
a1 <- ap[gene_id == "ENSG01.1" & chr == "chr1" & pos == 5000, slope / slope_se]
ok(isTRUE(all.equal(r1$alpha1, a1, tolerance = 1e-5)) && isTRUE(all.equal(r1$beta1, G1, tolerance = 1e-6)) &&
   isTRUE(all.equal(r1$beta2, G2, tolerance = 1e-6)), "alpha1, beta1, beta2 of the first feature")
cat(if (nfail == 0) "ALL TESTS PASSED\n" else sprintf("%d TESTS FAILED\n", nfail))
