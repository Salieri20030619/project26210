#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: install_infercnv_project_lib_v1.R <library_dir>")

lib <- normalizePath(args[[1L]], winslash = "/", mustWork = FALSE)
dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))

options(
  repos = c(CRAN = "https://cloud.r-project.org"),
  timeout = max(1800, getOption("timeout")),
  Ncpus = 8L
)

message("R_VERSION=", R.version.string)
message("TARGET_LIBRARY=", lib)
message("LIB_PATHS=", paste(.libPaths(), collapse = ";"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib, dependencies = TRUE)
}

installed <- FALSE
bioc_error <- NULL
tryCatch(
  {
    BiocManager::install("infercnv", lib = lib, ask = FALSE, update = FALSE, Ncpus = 8L)
    installed <- requireNamespace("infercnv", quietly = TRUE)
  },
  error = function(e) {
    bioc_error <<- conditionMessage(e)
    message("BIOCONDUCTOR_INSTALL_ERROR=", bioc_error)
  }
)

if (!installed) {
  if (!requireNamespace("remotes", quietly = TRUE)) {
    install.packages("remotes", lib = lib, dependencies = TRUE)
  }
  remotes::install_github(
    "broadinstitute/infercnv",
    lib = lib,
    dependencies = TRUE,
    upgrade = "never",
    build_vignettes = FALSE,
    Ncpus = 8L
  )
}

if (!requireNamespace("infercnv", quietly = TRUE)) stop("infercnv installation did not become loadable")
suppressPackageStartupMessages(library(infercnv))
message("INFERCNV_VERSION=", as.character(packageVersion("infercnv")))
message("INFERCNV_PATH=", find.package("infercnv"))

session <- file.path(lib, "INSTALL_SESSION_INFO_v1.txt")
writeLines(capture.output(sessionInfo()), session, useBytes = TRUE)
marker <- file.path(lib, "INSTALL_COMPLETE_v1")
part <- paste0(marker, ".part")
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), part, useBytes = TRUE)
if (!file.rename(part, marker)) stop("Atomic completion marker rename failed")
