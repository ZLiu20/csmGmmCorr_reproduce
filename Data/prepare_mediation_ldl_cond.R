here::i_am("Data/prepare_mediation_ldl_cond.R")
suppressPackageStartupMessages(library(data.table))

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")
  flush.console()
}
args <- commandArgs(trailingOnly = TRUE)
reference <- if (length(args)) args[1L] else Sys.getenv("PLINK_BFILE", "")
plink <- if (length(args) >= 2L) args[2L] else "plink"
bcftools <- if (length(args) >= 3L) args[3L] else "bcftools"
vcf_dir <- if (dir.exists(reference)) reference else NA_character_
if (is.na(vcf_dir) && !file.exists(paste0(reference, ".bim"))) {
  stop("Supply a PLINK binary prefix or a folder with 1000 Genomes VCFs")
}
rho2_proxy <- 0.8
root <- Sys.getenv("MEDLDL_OUT", here::here("Data"))
source_root <- Sys.getenv("MEDLDL_SRC", here::here("Data", "raw", "mediation_ldl"))
work_dir <- file.path(root, "prepared", "mediation_ldl", "conditional")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
input <- file.path(root, "mediation", "mediation_ldl_liver_noshare.tsv.gz")
if (!file.exists(input)) stop("Run Data/prepare_mediation_ldl.R first")
d <- fread(input)
genes <- unique(c(d$gene1, d$gene2, d$gene3))
say("features=%d, genes=%d", nrow(d), length(genes))

split_dir <- file.path(root, "prepared", "mediation_ldl", "liver_pairs")
lead_rows <- rbindlist(lapply(1:22, function(chromosome) {
  path <- file.path(split_dir, sprintf("chr%d.tsv", chromosome))
  if (!file.exists(path)) return(NULL)
  x <- fread(path, header = FALSE,
             col.names = c("gene", "variant", "pos", "ref", "alt", "maf", "beta", "se", "type"))
  x <- x[gene %chin% genes & is.finite(beta) & is.finite(se) & se > 0]
  x[, z := beta / se]
  x[, .SD[which.max(abs(z))], by = gene][, .(gene, lead = variant, lead_z = z)]
}))
if (!nrow(lead_rows)) stop("No lead eQTL variants were found")

lookup <- fread(cmd = paste("gzip -cd", shQuote(file.path(source_root, "GTEx_v8_lookup_table.txt.gz"))),
                select = c("variant_id", "variant_id_b37"))
needed <- unique(c(d$variant, lead_rows$lead))
lookup <- lookup[variant_id %chin% needed & !is.na(variant_id_b37) & nzchar(variant_id_b37)]
lookup[, c("chr", "pos", "ref", "alt") := tstrsplit(variant_id_b37, "_", fixed = TRUE, keep = 1:4)]
lookup[, `:=`(chr = as.integer(sub("^chr", "", chr, ignore.case = TRUE)), pos = as.integer(pos))]
lookup <- unique(lookup[!is.na(chr), .(variant = variant_id, chr, pos,
                                      ref = toupper(ref), alt = toupper(alt))])

if (!is.na(vcf_dir)) {
  panel_path <- file.path(vcf_dir, "integrated_call_samples_v3.20130502.ALL.panel")
  panel <- fread(panel_path)
  european <- panel[super_pop == "EUR", sample]
  sample_file <- file.path(work_dir, "eur_samples.txt")
  writeLines(european, sample_file)
  regions <- unique(lookup[, .(chr, pos)])[order(chr, pos)]
  parts <- character()
  for (chromosome in 1:22) {
    vcf <- file.path(vcf_dir, sprintf("ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz", chromosome))
    if (!file.exists(vcf)) stop("Missing reference VCF: ", vcf)
    region_file <- file.path(work_dir, sprintf("regions_chr%d.tsv", chromosome))
    fwrite(regions[chr == chromosome], region_file, sep = "\t", col.names = FALSE)
    part <- file.path(work_dir, sprintf("eur_chr%d.vcf.gz", chromosome))
    status <- system2(bcftools, c("view", "-R", region_file, "-S", sample_file,
                                  "-v", "snps", "-m2", "-M2", "-Oz", "-o", part, vcf))
    if (status != 0L) stop("bcftools view failed for chromosome ", chromosome)
    parts <- c(parts, part)
  }
  joined <- file.path(work_dir, "eur_all.vcf.gz")
  if (system2(bcftools, c("concat", "-Oz", "-o", joined, parts)) != 0L) stop("bcftools concat failed")
  reference <- file.path(work_dir, "g1000_eur_needed")
  if (system2(plink, c("--vcf", joined, "--double-id", "--set-missing-var-ids", "@:#[b37]$r,$a",
                       "--make-bed", "--out", reference)) != 0L) stop("PLINK VCF conversion failed")
}

bim <- fread(paste0(reference, ".bim"), header = FALSE,
             col.names = c("chr", "id", "cm", "pos", "a1", "a2"),
             colClasses = list(character = c("id", "a1", "a2")))
bim[, `:=`(a1 = toupper(a1), a2 = toupper(a2))]
matched <- merge(lookup, bim[, .(chr, pos, a1, a2)],
                 by = c("chr", "pos"), allow.cartesian = TRUE)
matched <- matched[(a1 == alt & a2 == ref) | (a1 == ref & a2 == alt)]
matched[, effect_sign := fifelse(a1 == alt, 1, -1)]
matched <- matched[!variant %chin% variant[duplicated(variant)]]
say("reference variants matched by position and alleles: %d", nrow(matched))

range_file <- file.path(work_dir, "range.tsv")
fwrite(matched[, .(chr, pos, pos, variant)], range_file, sep = "\t", col.names = FALSE)
subset_prefix <- file.path(work_dir, "reference_subset")
if (system2(plink, c("--bfile", reference, "--extract", "range", range_file,
                     "--make-bed", "--out", subset_prefix, "--allow-no-sex")) != 0L) {
  stop("PLINK reference subset failed")
}
subset_bim <- fread(paste0(subset_prefix, ".bim"), header = FALSE,
                    col.names = c("chr", "id", "cm", "pos", "a1", "a2"))
subset_bim[, `:=`(a1 = toupper(a1), a2 = toupper(a2))]
mapped_id <- merge(subset_bim, matched, by = c("chr", "pos", "a1", "a2"), all.x = TRUE)
if (anyNA(mapped_id$variant)) stop("Reference subset contains unmatched variants")
subset_bim[, id := mapped_id$variant[match(paste(chr, pos, a1, a2),
                                           paste(mapped_id$chr, mapped_id$pos,
                                                 mapped_id$a1, mapped_id$a2))]]
fwrite(subset_bim, paste0(subset_prefix, ".bim"), sep = "\t", col.names = FALSE)

lead_file <- file.path(work_dir, "lead_variants.txt")
writeLines(unique(lead_rows$lead), lead_file)
ld_prefix <- file.path(work_dir, "lead_ld")
if (system2(plink, c("--bfile", subset_prefix, "--r", "--ld-snp-list", lead_file,
                     "--ld-window-kb", "1500", "--ld-window", "999999",
                     "--ld-window-r2", "0", "--out", ld_prefix, "--allow-no-sex")) != 0L) {
  stop("PLINK correlation calculation failed")
}
ld <- fread(paste0(ld_prefix, ".ld"), select = c("SNP_A", "SNP_B", "R"))
sign_map <- setNames(matched$effect_sign, matched$variant)
ld[, R := R * sign_map[SNP_A] * sign_map[SNP_B]]
ld <- unique(rbind(ld[, .(lead = SNP_A, variant = SNP_B, r = R)],
                   ld[, .(lead = SNP_B, variant = SNP_A, r = R)]))
setkey(ld, lead, variant)

lead_by_gene <- setNames(lead_rows$lead, lead_rows$gene)
lead_z_by_gene <- setNames(lead_rows$lead_z, lead_rows$gene)
conditional_alpha <- function(variant, gene, z) {
  lead <- unname(lead_by_gene[[gene]])
  lead_z <- unname(lead_z_by_gene[[gene]])
  correlation <- ld[data.table(lead = lead, variant = variant), r]
  if (!length(correlation)) correlation <- NA_real_
  if (variant == lead) {
    value <- z; how <- "lead"
  } else if (!is.finite(correlation)) {
    value <- z; how <- "no_ld"
  } else if (correlation^2 >= rho2_proxy) {
    value <- z; how <- "proxy"
  } else {
    value <- (z - correlation * lead_z) / sqrt(1 - correlation^2)
    how <- "conditional"
  }
  list(value = value, how = how, r = correlation)
}

for (g in 1:3) {
  result <- Map(conditional_alpha, d$variant, d[[paste0("gene", g)]],
                d[[paste0("alpha", g)]])
  d[, paste0("alpha", g, "_marg") := get(paste0("alpha", g))]
  d[, paste0("alpha", g) := vapply(result, `[[`, numeric(1L), "value")]
  d[, paste0("alpha", g, "_how") := vapply(result, `[[`, character(1L), "how")]
  d[, paste0("r_lead", g) := vapply(result, `[[`, numeric(1L), "r")]
}

d <- d[complete.cases(d[, .(alpha1, alpha2, alpha3)])]
output <- file.path(root, "mediation", "mediation_ldl_liver_noshare_cond.tsv.gz")
fwrite(d, output, sep = "\t")
iqr_scale <- function(x) diff(quantile(x, c(.25, .75))) / 1.349
marginal <- as.matrix(d[, .(alpha1_marg, alpha2_marg, alpha3_marg)])
conditional <- as.matrix(d[, .(alpha1, alpha2, alpha3)])
print(data.table(
  alpha = 1:3,
  iqr_sd_marginal = round(apply(marginal, 2L, iqr_scale), 2L),
  iqr_sd_conditional = round(apply(conditional, 2L, iqr_scale), 2L),
  abs_gt4_marginal = round(colMeans(abs(marginal) > 4), 4L),
  abs_gt4_conditional = round(colMeans(abs(conditional) > 4), 4L)
))
say("DONE")
