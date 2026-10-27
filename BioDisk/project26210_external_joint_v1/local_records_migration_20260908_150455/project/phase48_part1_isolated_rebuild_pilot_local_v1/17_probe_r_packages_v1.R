#!/usr/bin/env Rscript

packages <- c(
  "Matrix", "data.table", "Seurat", "SeuratObject", "rhdf5", "hdf5r",
  "R.utils", "jsonlite", "digest"
)
cat("R_VERSION\t", R.version.string, "\n", sep = "")
for (package in packages) {
  present <- requireNamespace(package, quietly = TRUE)
  version <- if (present) as.character(utils::packageVersion(package)) else "NA"
  cat(package, present, version, sep = "\t")
  cat("\n")
}
