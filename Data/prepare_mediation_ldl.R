here::i_am("Data/prepare_mediation_ldl.R")
suppressPackageStartupMessages(library(data.table))

say <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), sprintf(...), "\n")
  flush.console()
}
args <- commandArgs(trailingOnly = TRUE)
source_root <- if (length(args)) args[1L] else here::here("Data", "raw", "mediation_ldl")
glgc_root <- if (length(args) >= 2L) args[2L] else here::here("Data", "raw", "glgc")
prepared_root <- Sys.getenv("MEDLDL_OUT", here::here("Data"))
work_dir <- file.path(prepared_root, "prepared", "mediation_ldl")
output_dir <- file.path(prepared_root, "mediation")
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

input <- function(name) file.path(source_root, name)
pair_file <- input("QTD000266.all.tsv.gz")
gtf_file <- input("gencode.v26.GRCh38.genes.gtf")
lookup_file <- input("GTEx_v8_lookup_table.txt.gz")
model_db <- input("mashr_liver.db")
covariance_file <- input("mashr_liver.txt.gz")
required <- c(pair_file, gtf_file, lookup_file, model_db, covariance_file)
missing <- required[!file.exists(required)]
if (length(missing)) stop("Missing mediation input: ", paste(missing, collapse = ", "))

find_glgc <- function(pattern) {
  hits <- list.files(glgc_root, pattern = pattern, recursive = TRUE,
                     full.names = TRUE, ignore.case = TRUE)
  if (!length(hits)) stop("GLGC LDL input not found for pattern: ", pattern)
  hits[1L]
}
ldl_files <- c(
  female = find_glgc("LDL.*EUR.*FEMALE.*singlevar.*[.]gz$"),
  male = find_glgc("LDL.*EUR.*MALE.*singlevar.*[.]gz$")
)

read_model_table <- function(table) {
  if (requireNamespace("RSQLite", quietly = TRUE)) {
    connection <- DBI::dbConnect(RSQLite::SQLite(), model_db)
    on.exit(DBI::dbDisconnect(connection))
    return(as.data.table(DBI::dbReadTable(connection, table)))
  }
  command <- sprintf("sqlite3 -header -csv %s %s",
                     shQuote(model_db), shQuote(sprintf("select * from %s;", table)))
  fread(cmd = command)
}

weights <- read_model_table("weights")
extra <- read_model_table("extra")
setnames(weights, tolower(names(weights)))
setnames(extra, tolower(names(extra)))
first_name <- function(candidates, available, label) {
  hit <- intersect(candidates, available)
  if (!length(hit)) stop("Cannot identify ", label, " column in MASHR model")
  hit[1L]
}
variant_col <- first_name(c("varid", "rsid", "variant", "id"), names(weights), "variant")
gene_col <- first_name(c("gene", "gene_id", "genename"), names(weights), "gene")
effect_col <- first_name(c("eff_allele", "effect_allele", "alt", "a1"), names(weights), "effect allele")
weight_col <- first_name(c("weight", "wgt", "value"), names(weights), "weight")
weights <- weights[, .(
  variant = as.character(get(variant_col)),
  gene = sub("[.][0-9]+$", "", as.character(get(gene_col))),
  effect = toupper(as.character(get(effect_col))),
  weight = as.numeric(get(weight_col))
)]
parts <- tstrsplit(weights$variant, "_", fixed = TRUE)
if (length(parts) < 4L) stop("MASHR variant identifiers must contain chromosome, position, REF and ALT")
weights[, `:=`(chr38 = sub("^chr", "", parts[[1L]], ignore.case = TRUE),
               pos38 = as.integer(parts[[2L]]), ref38 = toupper(parts[[3L]]),
               alt38 = toupper(parts[[4L]]))]
say("MASHR weights: %d rows across %d genes", nrow(weights), uniqueN(weights$gene))

lookup <- fread(cmd = paste("gzip -cd", shQuote(lookup_file)),
                select = c("variant_id", "variant_id_b37"))
lookup <- lookup[variant_id %chin% weights$variant & !is.na(variant_id_b37) & nzchar(variant_id_b37)]
lookup[, c("chr37", "pos37", "ref37", "alt37") :=
         tstrsplit(variant_id_b37, "_", fixed = TRUE, keep = 1:4)]
lookup[, key := paste(sub("^chr", "", chr37, ignore.case = TRUE),
                      as.integer(pos37), sep = ":")]

read_ldl <- function(path) {
  x <- fread(cmd = paste("gzip -cd", shQuote(path)),
             select = c("CHROM", "POS_b37", "REF", "ALT", "EFFECT_SIZE", "SE"),
             colClasses = list(character = c("CHROM", "REF", "ALT")))
  x[, key := paste(sub("^chr", "", CHROM, ignore.case = TRUE),
                   as.integer(POS_b37), sep = ":")]
  x[is.finite(EFFECT_SIZE) & is.finite(SE) & SE > 0,
    .(key, REF = toupper(REF), ALT = toupper(ALT), beta = EFFECT_SIZE, se = SE)]
}
female <- read_ldl(ldl_files[["female"]])
male <- read_ldl(ldl_files[["male"]])
ldl <- merge(female, male, by = c("key", "REF", "ALT"), all = TRUE,
             suffixes = c("_female", "_male"))
ldl[, `:=`(
  wf = fifelse(is.na(beta_female), 0, 1 / se_female^2),
  wm = fifelse(is.na(beta_male), 0, 1 / se_male^2)
)]
ldl[, z := (fifelse(is.na(beta_female), 0, beta_female * wf) +
              fifelse(is.na(beta_male), 0, beta_male * wm)) / sqrt(wf + wm)]
model_z <- merge(lookup[, .(variant = variant_id, key, ref37, alt37)],
                 ldl[, .(key, REF, ALT, z)], by = "key", allow.cartesian = TRUE)
model_z <- model_z[(REF == ref37 & ALT == alt37) | (REF == alt37 & ALT == ref37)]
model_z[, z_alt := fifelse(REF == ref37, z, -z)]
model_z <- model_z[!variant %chin% variant[duplicated(variant)]]
fwrite(model_z[, .(variant, z_alt)],
       file.path(work_dir, "ldl_model_variants.tsv.gz"), sep = "\t")
say("LDL statistics matched to %d model variants", nrow(model_z))

read_covariance <- function(path) {
  x <- fread(cmd = paste("gzip -cd", shQuote(path)))
  setnames(x, toupper(names(x)))
  needed <- c("GENE", "RSID1", "RSID2", "VALUE")
  if (length(setdiff(needed, names(x)))) stop("Unexpected MASHR covariance columns")
  x[, GENE := sub("[.][0-9]+$", "", GENE)]
  x
}

twas_path <- file.path(output_dir, "ldl_twas_liver.tsv.gz")
if (!file.exists(twas_path)) {
  covariance <- read_covariance(covariance_file)
  joined <- merge(weights, model_z[, .(variant, z = z_alt)], by = "variant")
  result <- rbindlist(lapply(split(joined, joined$gene), function(model) {
    gene <- model$gene[1L]
    variants <- model$variant
    cov <- covariance[GENE == gene & RSID1 %chin% variants & RSID2 %chin% variants]
    S <- matrix(0, length(variants), length(variants), dimnames = list(variants, variants))
    if (nrow(cov)) {
      S[cbind(match(cov$RSID1, variants), match(cov$RSID2, variants))] <- cov$VALUE
      S[cbind(match(cov$RSID2, variants), match(cov$RSID1, variants))] <- cov$VALUE
    }
    scale <- sqrt(pmax(diag(S), 0))
    variance <- drop(crossprod(model$weight, S %*% model$weight))
    if (!is.finite(variance) || variance <= 0) return(NULL)
    data.table(gene = gene, n_model = nrow(model), n_used = sum(is.finite(model$z)),
               z_twas = sum(model$weight * scale * model$z) / sqrt(variance))
  }), fill = TRUE)
  if (!nrow(result)) stop("No S-PrediXcan statistics could be computed")
  gene_name_col <- intersect(c("gene_name", "genename", "name"), names(extra))[1L]
  extra_gene_col <- intersect(c("gene", "gene_id"), names(extra))[1L]
  if (!is.na(gene_name_col) && !is.na(extra_gene_col)) {
    annotation <- extra[, .(gene = sub("[.][0-9]+$", "", as.character(get(extra_gene_col))),
                            gene_name = as.character(get(gene_name_col)))]
    result <- merge(result, unique(annotation), by = "gene", all.x = TRUE)
  } else result[, gene_name := gene]
  fwrite(result, twas_path, sep = "\t")
}
twas <- fread(twas_path)

gtf <- fread(cmd = paste(if (grepl("[.]gz$", gtf_file)) "gzip -cd" else "cat",
                         shQuote(gtf_file), "| grep -v '^#'"),
             sep = "\t", header = FALSE, select = c(1, 3, 4, 5, 7, 9),
             col.names = c("chr", "feature", "start", "end", "strand", "attributes"))
gtf <- gtf[feature == "gene"]
gtf[, gene := sub('.*gene_id "([^"]+)".*', "\\1", attributes)]
gtf[, gene := sub("[.][0-9]+$", "", gene)]
gtf[, tss := fifelse(strand == "-", end, start)]
gtf <- unique(gtf[gene %chin% twas$gene, .(gene, tss)])
say("Liver prediction genes with annotated TSS: %d", nrow(gtf))

split_dir <- file.path(work_dir, "liver_pairs")
done_marker <- file.path(split_dir, "DONE")
if (!file.exists(done_marker)) {
  dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
  wanted <- file.path(work_dir, "model_genes.tsv")
  fwrite(gtf[, .(gene)], wanted, sep = "\t", col.names = FALSE)
  header <- names(fread(cmd = paste("gzip -cd", shQuote(pair_file), "| head -1"), nrows = 0L))
  aliases <- list(
    gene = c("gene_id", "molecular_trait_id"), chromosome = c("chromosome", "chr"),
    position = c("position", "pos"), ref = c("ref", "reference_allele"),
    alt = c("alt", "alternate_allele"), variant = c("variant", "variant_id"),
    maf = c("maf", "minor_allele_frequency"), beta = c("beta", "slope"),
    se = c("se", "slope_se"), type = c("type", "variant_type")
  )
  column <- vapply(aliases, function(x) {
    i <- match(tolower(x), tolower(header), nomatch = 0L)
    i <- i[i > 0L]
    if (!length(i)) stop("Missing eQTL Catalogue column: ", paste(x, collapse = "/"))
    i[1L]
  }, integer(1L))
  awk <- file.path(work_dir, "split_eqtl.awk")
  writeLines(c(
    "BEGIN { FS=OFS=\"\\t\"; while ((getline g < genes) > 0) keep[g]=1 }",
    sprintf("NR>1 { gene=$%d; sub(/\\.[0-9]+$/,\"\",gene); if (!(gene in keep)) next", column[["gene"]]),
    sprintf("chr=$%d; sub(/^chr/,\"\",chr); if (chr !~ /^[0-9]+$/) next", column[["chromosome"]]),
    sprintf("if ($%d+0 < 0.05 || $%d != \"SNP\") next", column[["maf"]], column[["type"]]),
    sprintf("print gene,$%d,$%d,$%d,$%d,$%d,$%d,$%d,$%d > (out \"/chr\" chr \".tsv\") }",
            column[["variant"]], column[["position"]], column[["ref"]], column[["alt"]],
            column[["maf"]], column[["beta"]], column[["se"]], column[["type"]])
  ), awk)
  status <- system2("bash", c("-c", shQuote(sprintf("gzip -cd %s | awk -v genes=%s -v out=%s -f %s",
                                                      shQuote(pair_file), shQuote(wanted),
                                                      shQuote(split_dir), shQuote(awk)))))
  if (status != 0L) stop("Failed to split the eQTL Catalogue file")
  writeLines("DONE", done_marker)
}

twas_z <- setNames(twas$z_twas, twas$gene)
twas_name <- setNames(twas$gene_name, twas$gene)
model_sets <- split(weights[, .(variant)], weights$gene)
shares_variant <- function(a, b) length(intersect(model_sets[[a]]$variant, model_sets[[b]]$variant))
choose_genes <- function(genes, rule, number = 3L) {
  if (length(genes) < number) return(NULL)
  if (rule == "nearest") return(genes[seq_len(number)])
  selected <- genes[1L]
  for (gene in genes[-1L]) {
    if (all(vapply(selected, function(old) shares_variant(old, gene) == 0L, logical(1L)))) {
      selected <- c(selected, gene)
    }
    if (length(selected) == number) return(selected)
  }
  NULL
}

outputs <- list(noshare = list(), nearest = list())
for (chromosome in 1:22) {
  path <- file.path(split_dir, sprintf("chr%d.tsv", chromosome))
  if (!file.exists(path)) next
  eqtl <- fread(path, header = FALSE,
                col.names = c("gene", "variant", "pos", "ref", "alt", "maf", "beta", "se", "type"))
  eqtl <- eqtl[is.finite(beta) & is.finite(se) & se > 0]
  eqtl[, z := beta / se]
  eqtl <- merge(eqtl, gtf, by = "gene")
  variants <- unique(eqtl[, .(variant, pos)])
  variants[, window := pos %/% 10000L]
  variants[, midpoint_distance := abs(pos - (window * 10000L + 5000L))]
  setorder(variants, window, midpoint_distance, variant)
  variants <- variants[, .SD[1L], by = window]
  candidate <- eqtl[variant %chin% variants$variant]
  candidate[, tss_distance := abs(tss - pos)]
  setorder(candidate, variant, tss_distance, gene)
  by_variant <- split(candidate, candidate$variant)
  for (rule in names(outputs)) {
    rows <- lapply(by_variant, function(x) {
      genes <- unique(x$gene)
      selected <- choose_genes(genes, rule)
      if (is.null(selected)) return(NULL)
      row <- list(variant = x$variant[1L], chr = chromosome, pos_b38 = x$pos[1L],
                  shared = sum(combn(selected, 2L, function(g) shares_variant(g[1L], g[2L]))))
      for (g in seq_along(selected)) {
        hit <- x[gene == selected[g]][1L]
        row[[paste0("gene", g)]] <- selected[g]
        row[[paste0("gene", g, "_name")]] <- twas_name[[selected[g]]]
        row[[paste0("tss", g)]] <- hit$tss
        row[[paste0("alpha", g)]] <- hit$z
        row[[paste0("beta", g)]] <- twas_z[[selected[g]]]
      }
      as.data.table(row)
    })
    outputs[[rule]][[length(outputs[[rule]]) + 1L]] <- rbindlist(rows, fill = TRUE)
  }
  say("chromosome %d: %d candidate variants", chromosome, nrow(candidate))
}

stat_columns <- as.vector(rbind(paste0("alpha", 1:3), paste0("beta", 1:3)))
for (rule in names(outputs)) {
  result <- rbindlist(outputs[[rule]], fill = TRUE)
  setorder(result, chr, pos_b38, variant)
  result <- result[complete.cases(result[, ..stat_columns])]
  fwrite(result, file.path(output_dir, sprintf("mediation_ldl_liver_%s.tsv.gz", rule)), sep = "\t")
  say("%s: %d features and %d gene triples", rule, nrow(result),
      uniqueN(result[, paste(gene1, gene2, gene3, sep = ":")]))
}
say("DONE")
