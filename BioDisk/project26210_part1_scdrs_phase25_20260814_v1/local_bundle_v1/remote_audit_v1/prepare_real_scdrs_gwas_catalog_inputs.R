#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(jsonlite)
})

out_dir <- Sys.getenv("OUT_DIR", "results/scDRS_real_GWAS_20260713")
input_rds <- Sys.getenv("INPUT_RDS", "object/11_malignant_epithelial_cNMF_inferCNV_copyKAT_integrated.rds")
assoc_dir <- file.path(out_dir, "downloads", "associations")
selected_file <- file.path(out_dir, "downloads", "RCC_susceptibility_GWAS_selected.tsv")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "inputs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)

message("[scDRS-input] Reading malignant epithelial Seurat object: ", input_rds)
obj <- readRDS(input_rds)
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
if (nrow(counts) == 0) counts <- GetAssayData(obj, assay = "RNA", slot = "counts")
counts <- as(counts, "dgCMatrix")

message("[scDRS-input] Parsing GWAS Catalog associations")
assoc_files <- Sys.glob(file.path(assoc_dir, "*.associations.json"))
if (!length(assoc_files)) stop("No association JSON files found in: ", assoc_dir)

rows <- list()
k <- 1L
for (f in assoc_files) {
  acc <- sub("\\.associations\\.json$", "", basename(f))
  x <- fromJSON(f, simplifyVector = FALSE)
  associations <- x[["_embedded"]][["associations"]]
  if (is.null(associations) || !length(associations)) next
  for (a in associations) {
    p <- suppressWarnings(as.numeric(a[["pvalue"]]))
    if (!is.finite(p) || p <= 0) next
    loci <- a[["loci"]]
    if (is.null(loci) || !length(loci)) next
    for (loc in loci) {
      genes <- loc[["authorReportedGenes"]]
      if (is.null(genes) || !length(genes)) next
      for (g in genes) {
        gene <- g[["geneName"]]
        if (is.null(gene) || is.na(gene) || !nzchar(gene)) next
        rows[[k]] <- data.table(
          study_accession = acc,
          gene = gene,
          pvalue = p,
          log10p = -log10(p)
        )
        k <- k + 1L
      }
    }
  }
}
gene_assoc <- unique(rbindlist(rows, fill = TRUE))
gene_assoc <- gene_assoc[gene %in% rownames(counts)]
if (nrow(gene_assoc) < 10) {
  stop("Too few GWAS Catalog genes present in expression matrix: ", nrow(gene_assoc))
}

gene_summary <- gene_assoc[
  ,
  .(
    n_assoc = .N,
    n_studies = uniqueN(study_accession),
    min_pvalue = min(pvalue, na.rm = TRUE),
    max_log10p = max(log10p, na.rm = TRUE),
    studies = paste(sort(unique(study_accession)), collapse = ",")
  ),
  by = gene
]
gene_summary[, z_weight := pmin(qnorm(min_pvalue / 2, lower.tail = FALSE), 10)]
gene_summary <- gene_summary[order(-z_weight, gene)]

gs_trait <- "RCC_GWAS_Catalog_loci"
gs_string <- paste0(gene_summary$gene, ":", signif(gene_summary$z_weight, 4), collapse = ",")
gs <- data.table(TRAIT = gs_trait, GENESET = gs_string)

fwrite(gene_assoc, file.path(out_dir, "tables", "RCC_GWAS_Catalog_loci_gene_associations.tsv"), sep = "\t")
fwrite(gene_summary, file.path(out_dir, "tables", "RCC_GWAS_Catalog_loci_gene_summary.tsv"), sep = "\t")
fwrite(gs, file.path(out_dir, "inputs", "RCC_GWAS_Catalog_loci.gs"), sep = "\t")

message("[scDRS-input] Writing sparse counts MatrixMarket as cells x genes")
writeMM(t(counts), file.path(out_dir, "inputs", "malignant_epithelial_counts_cells_by_genes.mtx"))
fwrite(data.table(gene = rownames(counts)), file.path(out_dir, "inputs", "malignant_epithelial_genes.tsv"), sep = "\t")

meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
n_counts <- Matrix::colSums(counts)
n_genes <- Matrix::colSums(counts > 0)
cov <- data.table(
  index = colnames(counts),
  const = 1,
  n_genes = as.numeric(n_genes),
  log_n_counts = log1p(as.numeric(n_counts))
)
meta <- merge(data.table(cell = colnames(counts)), meta, by = "cell", all.x = TRUE, sort = FALSE)
fwrite(meta, file.path(out_dir, "inputs", "malignant_epithelial_obs.tsv"), sep = "\t")
fwrite(cov, file.path(out_dir, "inputs", "malignant_epithelial_scdrs.cov"), sep = "\t")

selected <- if (file.exists(selected_file)) fread(selected_file) else data.table()
manifest <- data.table(
  item = c(
    "input_rds",
    "n_cells",
    "n_genes",
    "n_gwas_studies_selected",
    "n_gwas_genes_in_expression",
    "gs_file",
    "h5ad_expected",
    "cov_file",
    "method_note"
  ),
  value = c(
    input_rds,
    ncol(counts),
    nrow(counts),
    if (nrow(selected)) uniqueN(selected[["STUDY ACCESSION"]]) else NA_integer_,
    nrow(gene_summary),
    file.path(out_dir, "inputs", "RCC_GWAS_Catalog_loci.gs"),
    file.path(out_dir, "inputs", "malignant_epithelial_counts.h5ad"),
    file.path(out_dir, "inputs", "malignant_epithelial_scdrs.cov"),
    "GWAS Catalog significant author-reported genes from RCC/kidney cancer susceptibility studies; weights are two-sided z from minimum association p-value. This is loci-derived scDRS, not MAGMA whole-summary-stat scDRS."
  )
)
fwrite(manifest, file.path(out_dir, "inputs", "scDRS_real_GWAS_input_manifest.tsv"), sep = "\t")
message("[scDRS-input] Genes in RCC GWAS loci gene set: ", nrow(gene_summary))
message("[scDRS-input] Done.")
