#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mistyR)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: smoke_test_misty_deldir_fallback_v1.R MEASUREMENTS POSITIONS")
gate <- "/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh"
if (system2(gate) != 0L) stop("Memory gate failed")

measurements_raw <- read_tsv(args[[1]], show_col_types = FALSE)
positions_raw <- read_tsv(args[[2]], show_col_types = FALSE)
if (!identical(measurements_raw$barcode, positions_raw$barcode)) stop("Barcode mismatch")
measurements <- as.data.frame(measurements_raw[, -1, drop = FALSE])
rownames(measurements) <- measurements_raw$barcode
positions <- as.data.frame(positions_raw[, c("x_um", "y_um")])
rownames(positions) <- positions_raw$barcode

epsilon <- 1e-3
indices <- seq_len(nrow(positions))
jittered <- positions
jittered$x_um <- positions$x_um + epsilon * sin(indices * sqrt(2))
jittered$y_um <- positions$y_um + epsilon * cos(indices * sqrt(3))
max_shift <- max(abs(as.matrix(jittered) - as.matrix(positions)))
if (max_shift > epsilon * (1 + 1e-12)) stop("Unexpected jitter magnitude")

views <- create_initial_view(measurements, unique.id = "deldir_fallback_smoke_v1")
views <- add_juxtaview(views, jittered, neighbor.thr = 110, cached = FALSE, verbose = FALSE)
views <- add_paraview(views, jittered, l = 200, zoi = 110, family = "gaussian", cached = FALSE, verbose = FALSE)
view_names <- names(views)
if (!any(grepl("juxta", view_names)) || !any(grepl("para", view_names))) {
  stop("Expected juxta/para views were not constructed")
}
cat(sprintf(
  "MISTY_DELDIR_FALLBACK_SMOKE_PASS spots=%d features=%d epsilon_um=%.6g max_axis_shift_um=%.6g views=%s\n",
  nrow(measurements), ncol(measurements), epsilon, max_shift, paste(view_names, collapse = ",")
))
