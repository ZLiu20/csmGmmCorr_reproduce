# Cis restriction of the lung SCC mediation data (Data/prepare_mediation_scc.R must be run first).
# A SNP-gene pair is kept when the gene is on the same chromosome as the SNP and its transcription
# start site (TSS) is within 1 Mb of the SNP (the GTEx cis window). Gene positions (GRCh37) come from
# ensembl_refgene_hg19_20180109.rda, the annotation used by Sun, McCaw and Lin (2025, Fig4 code).
# Transcripts are collapsed to genes by HGNC name; a gene annotated on more than one chromosome is
# dropped, and a gene with several TSSs keeps the TSS closest to the SNP.
# Output: Data/mediation/mediation_scc_cis.tsv.gz (same columns plus tss_dist_kb),
#         Data/mediation/mediation_scc_cis_counts.tsv
#   Rscript Data/prepare_mediation_scc_cis.R        (from the repository root)
here::i_am("Data/prepare_mediation_scc_cis.R")
suppressMessages(library(data.table))
WINDOW <- 1e6
dd <- here::here("Data", "mediation")
ann_file <- Sys.getenv("REFGENE_RDA", here::here("Data", "ensembl_refgene_hg19_20180109.rda"))
stopifnot(file.exists(ann_file))
e <- new.env(); load(ann_file, envir = e)
ann <- as.data.table(get(ls(e)[1], envir = e))
need <- c("HGNC_name", "Chr", "Strand", "txStart", "txEnd")
stopifnot(all(need %in% names(ann)))
ann <- ann[!is.na(HGNC_name) & HGNC_name != "", .(gene = HGNC_name, chr = sub("^chr", "", as.character(Chr)),
                                                  strand = as.character(Strand), txStart = as.numeric(txStart), txEnd = as.numeric(txEnd))]
ann <- ann[chr %in% as.character(1:22)]
ann[, chr := as.integer(chr)]
ann[, tss := ifelse(strand %in% c("-", "-1"), txEnd, txStart)]
multi <- ann[, uniqueN(chr), by = gene][V1 > 1, gene]
ann <- unique(ann[!gene %in% multi, .(gene, chr, tss)])
cat(sprintf("annotation: %d genes on autosomes (%d genes on >1 chromosome dropped)\n", uniqueN(ann$gene), length(multi)))

med <- fread(file.path(dd, "mediation_scc.tsv.gz"))
cat(sprintf("input: %d SNP-gene pairs, %d genes; genes without annotation: %d\n",
            nrow(med), uniqueN(med$gene), length(setdiff(unique(med$gene), ann$gene))))
x <- merge(med, ann, by = c("gene", "chr"), allow.cartesian = TRUE)
x[, tss_dist_kb := round((tss - pos_b37) / 1e3, 1)]
x <- x[abs(tss - pos_b37) <= WINDOW]
x <- x[order(abs(tss_dist_kb))]
x <- x[!duplicated(x, by = c("snp", "gene"))]
x[, tss := NULL]
setcolorder(x, c("snp", "locus_gene", "chr", "pos_b37", "gene", "tss_dist_kb",
                 "alpha_lung", "beta_lung", "alpha_blood", "beta_blood", "alpha_liver", "beta_liver"))
x <- x[order(chr, pos_b37, tss_dist_kb)]
fwrite(x, file.path(dd, "mediation_scc_cis.tsv.gz"), sep = "\t")
rank <- fread(file.path(dd, "mediation_scc_snp_rank.tsv"))[, .(snp, rank_ILCCO, zILCCO)]
dup <- x[, .N, by = gene][N > 1, gene]
cnt <- merge(x[, .(cis_pairs = .N, genes_shared_with_other_snp = sum(gene %in% dup)), by = .(snp, locus_gene, chr, pos_b37)], rank, by = "snp")
cnt <- cnt[order(rank_ILCCO)]
fwrite(cnt, file.path(dd, "mediation_scc_cis_counts.tsv"), sep = "\t")
options(width = 200); print(cnt)
cat(sprintf("cis: %d SNP-gene pairs, %d SNPs, %d genes (%d genes paired with >1 SNP)\n",
            nrow(x), uniqueN(x$snp), uniqueN(x$gene), length(dup)))
top3 <- head(rank[order(rank_ILCCO)]$snp, 3)
cat(sprintf("cis, top 3 SNPs (%s): %d pairs\n", paste(top3, collapse = ", "), nrow(x[snp %in% top3])))
