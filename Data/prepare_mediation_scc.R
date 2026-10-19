# Lung squamous cell carcinoma (SCC) mediation data for the correlated csmGmm.
# Feature = SNP-gene pair. Pathways = three tissues (lung, whole blood, liver), so G = 3 and
# K = 6 with coordinates (alpha_lung, beta_lung, alpha_blood, beta_blood, alpha_liver, beta_liver).
#   alpha = GTEx v8 eQTL statistic of the SNP for the expression of the gene in the tissue
#   beta  = S-PrediXcan statistic of the predicted expression of the gene in the tissue for SCC
# Inputs (data provided by R. Sun):
#   $ZILI/data/<rs>_<GENE>_<chr>_<pos>.txt          lung eQTL (15 SNPs)
#   $ZILI/data/scc_lung_addchr1.csv                  lung S-PrediXcan
#   $ZILI/zili_new_mediation_data.zip               whole-blood and liver eQTL and S-PrediXcan
#   $ZILI/data/lc_overall_three.txt                 lung cancer GWAS (only to rank the SNPs)
# Output: Data/mediation/mediation_scc.tsv.gz, Data/mediation/mediation_scc_counts.tsv,
#         Data/mediation/mediation_scc_snp_rank.tsv
#   Rscript Data/prepare_mediation_scc.R           (from the repository root)
here::i_am("Data/prepare_mediation_scc.R")
suppressMessages(library(data.table))
ZILI <- Sys.getenv("ZILI_DIR", "/rsrch3/home/biostatistics/UK_Biobank/zili")
od <- here::here("Data", "mediation"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
raw <- file.path(od, "scc_raw"); dir.create(raw, showWarnings = FALSE)
zipdir <- file.path(raw, "zili_new_mediation_data")
if (!dir.exists(zipdir)) utils::unzip(file.path(ZILI, "zili_new_mediation_data.zip"), exdir = raw)
stopifnot(dir.exists(zipdir))

# the 15 SNPs, from the liver eQTL file names: <rs>_<GENE>_<chr>_<pos>_Liver.txt
stems <- sub("_Liver\\.txt$", "", list.files(zipdir, pattern = "^rs.*_Liver\\.txt$"))
stopifnot(length(stems) == 15)
parts <- tstrsplit(stems, "_")
snps <- data.table(stem = stems, snp = parts[[1]], locus_gene = parts[[2]],
                   chr = as.integer(parts[[3]]), pos_b37 = as.integer(parts[[4]]))

read_eqtl <- function(f) {
  d <- fread(f, select = c("Gene", "testStat"))
  d <- d[!is.na(testStat) & !(duplicated(Gene) | duplicated(Gene, fromLast = TRUE))]
  setnames(d, c("gene", "z"))[]
}
read_twas <- function(f) {
  d <- fread(f, select = c("gene_name", "zscore"))
  d <- d[!is.na(zscore) & !(duplicated(gene_name) | duplicated(gene_name, fromLast = TRUE))]
  setnames(d, c("gene", "z"))[]
}
beta <- list(lung  = read_twas(file.path(ZILI, "data", "scc_lung_addchr1.csv")),
             blood = read_twas(file.path(zipdir, "scc__PM__Whole_Blood.csv")),
             liver = read_twas(file.path(zipdir, "scc__PM__Liver.csv")))
cat(sprintf("S-PrediXcan genes: lung %d, blood %d, liver %d\n", nrow(beta$lung), nrow(beta$blood), nrow(beta$liver)))

out <- list(); counts <- list()
for (i in seq_len(nrow(snps))) {
  s <- snps[i]
  alpha <- list(lung  = read_eqtl(file.path(ZILI, "data", paste0(s$stem, ".txt"))),
                blood = read_eqtl(file.path(zipdir, paste0(s$stem, "_Whole_Blood.txt"))),
                liver = read_eqtl(file.path(zipdir, paste0(s$stem, "_Liver.txt"))))
  d <- NULL
  for (tis in c("lung", "blood", "liver")) {
    x <- merge(alpha[[tis]], beta[[tis]], by = "gene", suffixes = c("_a", "_b"))
    setnames(x, c("z_a", "z_b"), paste0(c("alpha_", "beta_"), tis))
    d <- if (is.null(d)) x else merge(d, x, by = "gene")
  }
  d <- cbind(s[, .(snp, locus_gene, chr, pos_b37)], d)
  counts[[i]] <- data.table(snp = s$snp, locus_gene = s$locus_gene,
                            eqtl_lung = nrow(alpha$lung), eqtl_blood = nrow(alpha$blood), eqtl_liver = nrow(alpha$liver),
                            pairs_all_three_tissues = nrow(d))
  out[[i]] <- d
}
med <- rbindlist(out)
setcolorder(med, c("snp", "locus_gene", "chr", "pos_b37", "gene",
                   "alpha_lung", "beta_lung", "alpha_blood", "beta_blood", "alpha_liver", "beta_liver"))
fwrite(med, file.path(od, "mediation_scc.tsv.gz"), sep = "\t")
fwrite(rbindlist(counts), file.path(od, "mediation_scc_counts.tsv"), sep = "\t")
print(rbindlist(counts))
cat(sprintf("Total: %d SNP-gene pairs, %d SNPs, %d genes\n", nrow(med), uniqueN(med$snp), uniqueN(med$gene)))

# Rank the SNPs by their lung cancer association in ILCCO (external to the alpha and beta
# statistics), so that a reduced SNP set can be prespecified without looking at the analysis.
gw <- file.path(ZILI, "data", "lc_overall_three.txt")
hdr <- names(fread(gw, nrows = 0))
g <- fread(cmd = paste("grep -w -E", shQuote(paste(snps$snp, collapse = "|")), shQuote(gw)), header = FALSE)
setnames(g, hdr)
g <- g[RS %in% snps$snp, .(snp = RS, zILCCO, zUKB, zMVP)][!duplicated(snp)]
rank <- merge(snps[, .(snp, locus_gene, chr, pos_b37)], g, by = "snp", all.x = TRUE)[order(-abs(zILCCO))]
rank[, rank_ILCCO := seq_len(.N)]
fwrite(rank, file.path(od, "mediation_scc_snp_rank.tsv"), sep = "\t")
print(rank)
