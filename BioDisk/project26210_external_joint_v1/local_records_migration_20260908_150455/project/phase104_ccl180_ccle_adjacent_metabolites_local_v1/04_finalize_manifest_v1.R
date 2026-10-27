#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(digest)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: 04_finalize_manifest_v1.R <phase104_output_dir>")
root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
manifest_path <- file.path(root, "FILE_MANIFEST_SHA256_v1.tsv")
files <- list.files(root, recursive = TRUE, full.names = TRUE)
files <- files[normalizePath(files, winslash = "/") != normalizePath(manifest_path, winslash = "/", mustWork = FALSE)]
prefix <- paste0(root, "/")
manifest <- data.table(
  relative_path = substring(normalizePath(files, winslash = "/"), nchar(prefix) + 1L),
  bytes = file.info(files)$size,
  sha256 = vapply(files, digest, character(1), algo = "sha256", file = TRUE)
)
setorder(manifest, relative_path)
fwrite(manifest, manifest_path, sep = "\t", quote = FALSE)
cat("PASS_FINAL_MANIFEST\n")
cat("files=", nrow(manifest), "\n", sep = "")
