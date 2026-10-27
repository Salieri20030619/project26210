args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Usage: Rscript 03_audit_emtab3218_expression_v1.R <rma_matrix> <annotation> <sample_manifest> <signature_genes> <output_dir>")
}

rma_path <- args[[1L]]
annotation_path <- args[[2L]]
sample_manifest_path <- args[[3L]]
signature_path <- args[[4L]]
output_dir <- args[[5L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(frame, path) {
  write.table(frame, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

expr <- as.matrix(read.delim(
  rma_path, header = TRUE, row.names = 1L, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "", comment.char = ""
))
storage.mode(expr) <- "double"
annotation <- read.delim(
  annotation_path, header = TRUE, check.names = FALSE,
  stringsAsFactors = FALSE, quote = "", comment.char = "", fill = TRUE
)
samples <- read.delim(sample_manifest_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
signatures <- read.delim(signature_path, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
signatures$weight <- signatures$raw_weight

required_annotation <- c("Probeset", "Symbol")
required_samples <- c("sample_id", "array_file")
required_signature <- c("signature_id", "gene", "weight")
if (length(setdiff(required_annotation, names(annotation)))) stop("Annotation columns missing")
if (length(setdiff(required_samples, names(samples)))) stop("Sample manifest columns missing")
if (length(setdiff(required_signature, names(signatures)))) stop("Signature columns missing")

samples$matrix_column <- paste0("ENTREZG-", sub("\\.cel$", "", samples$array_file, ignore.case = TRUE))
samples$matrix_column_present <- samples$matrix_column %in% colnames(expr)
samples$matrix_column_index <- match(samples$matrix_column, colnames(expr))

f1 <- signatures[signatures$signature_id == "f1_top100_weighted", , drop = FALSE]
f1$weight <- as.numeric(f1$weight)
symbol_to_probeset <- split(annotation$Probeset, annotation$Symbol)
f1$probesets <- vapply(f1$gene, function(gene) {
  hits <- intersect(symbol_to_probeset[[gene]], rownames(expr))
  paste(hits, collapse = ";")
}, character(1L))
f1$present <- nzchar(f1$probesets)

gates <- data.frame(
  gate = c(
    "rma_loci", "rma_samples", "unique_rma_rows", "unique_rma_columns",
    "sample_manifest_rows", "sample_columns_mapped", "unique_sample_mapping",
    "f1_signature_genes", "f1_genes_present", "tymp_present", "finite_expression"
  ),
  observed = c(
    nrow(expr), ncol(expr), !anyDuplicated(rownames(expr)), !anyDuplicated(colnames(expr)),
    nrow(samples), sum(samples$matrix_column_present), !anyDuplicated(samples$matrix_column_index),
    nrow(f1), sum(f1$present), "TYMP" %in% annotation$Symbol,
    all(is.finite(expr))
  ),
  expected_or_threshold = c(18562, 114, TRUE, TRUE, 114, 114, TRUE, 100, 70, TRUE, TRUE),
  stringsAsFactors = FALSE
)
gates$status <- c(
  ifelse(nrow(expr) >= 18500L, "PASS", "FAIL"),
  ifelse(ncol(expr) == 114L, "PASS", "FAIL"),
  ifelse(!anyDuplicated(rownames(expr)), "PASS", "FAIL"),
  ifelse(!anyDuplicated(colnames(expr)), "PASS", "FAIL"),
  ifelse(nrow(samples) == 114L, "PASS", "FAIL"),
  ifelse(sum(samples$matrix_column_present) == 114L, "PASS", "FAIL"),
  ifelse(!anyDuplicated(samples$matrix_column_index), "PASS", "FAIL"),
  ifelse(nrow(f1) == 100L, "PASS", "FAIL"),
  ifelse(sum(f1$present) >= 70L, "PASS", "FAIL"),
  ifelse("TYMP" %in% annotation$Symbol, "PASS", "FAIL"),
  ifelse(all(is.finite(expr)), "PASS", "FAIL")
)

write_tsv(gates, file.path(output_dir, "expression_gate_v1.tsv"))
write_tsv(samples, file.path(output_dir, "expression_sample_mapping_v1.tsv"))
write_tsv(f1, file.path(output_dir, "f1_gene_coverage_v1.tsv"))
write_tsv(
  data.frame(
    metric = c("rma_min", "rma_median", "rma_max"),
    value = c(min(expr), median(expr), max(expr))
  ),
  file.path(output_dir, "expression_distribution_v1.tsv")
)

print(gates, row.names = FALSE)
if (any(gates$status == "FAIL")) quit(status = 2L)
cat("PASS_E_MTAB_3218_EXPRESSION_GATE\n")
