args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 00_preflight_r_packages_v1.R output.tsv")
packages <- c("CellChat", "nichenetr", "Seurat", "SeuratObject")
rows <- data.frame(
  package = packages,
  installed = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(packages, function(x) if (requireNamespace(x, quietly = TRUE)) as.character(packageVersion(x)) else "NOT_INSTALLED", character(1)),
  stringsAsFactors = FALSE
)
out <- args[[1]]
if (file.exists(out)) stop("Refusing to overwrite output")
part <- paste0(out, ".part")
write.table(rows, part, sep = "\t", quote = FALSE, row.names = FALSE)
if (!file.rename(part, out)) stop("Atomic rename failed")
