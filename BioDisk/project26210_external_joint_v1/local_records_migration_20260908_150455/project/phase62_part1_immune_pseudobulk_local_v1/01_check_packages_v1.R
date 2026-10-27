#!/usr/bin/env Rscript
packages <- c('Seurat', 'data.table', 'Matrix', 'edgeR', 'limma', 'fgsea', 'msigdbr')
result <- data.frame(
  package = packages,
  available = vapply(packages, requireNamespace, logical(1), quietly = TRUE),
  version = vapply(packages, function(x) if (requireNamespace(x, quietly = TRUE)) as.character(packageVersion(x)) else NA_character_, character(1))
)
print(result, row.names = FALSE)
