#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(mistyR)
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(future)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: run_part4_misty_sample_v1.R SAMPLE INPUT_DIR OUTPUT_DIR SEED")
}

sample_id <- args[[1]]
input_dir <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
seed <- as.integer(args[[4]])
if (!grepl("^[A-Za-z0-9_]+$", sample_id) || is.na(seed)) {
  stop("Invalid sample or seed")
}

project_root <- "/mnt/BioDisk/ludingyu/project26210"
memory_gate <- file.path(project_root, "check_memory_gate.sh")
measurement_path <- file.path(input_dir, paste0(sample_id, "_misty_measurements_v1.tsv"))
position_path <- file.path(input_dir, paste0(sample_id, "_misty_positions_v1.tsv"))
feature_manifest_path <- file.path(input_dir, "misty_feature_manifest_v1.tsv")
raw_dir <- file.path(output_dir, "raw_misty")

md5sum_file <- function(path) {
  unname(tools::md5sum(path))
}

write_status <- function(status, message = NULL) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  payload <- list(
    status = status,
    sample = sample_id,
    message = message,
    seed = seed,
    updated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  write_json(payload, file.path(output_dir, "status.json"), pretty = TRUE, auto_unbox = TRUE)
}

run_memory_gate <- function(stage) {
  message("MEMORY_GATE stage=", stage)
  output <- system2(memory_gate, stdout = TRUE, stderr = TRUE)
  status <- attr(output, "status")
  if (is.null(status)) status <- 0L
  cat(paste(output, collapse = "\n"), "\n")
  if (status != 0L) stop("Memory gate failed at ", stage)
}

safe_label <- function(values, manifest) {
  labels <- setNames(manifest$source_label, manifest$feature_id)
  result <- unname(labels[values])
  result[is.na(result)] <- values[is.na(result)]
  result
}

main <- function() {
  if (file.exists(file.path(output_dir, "status.json"))) {
    prior <- tryCatch(read_json(file.path(output_dir, "status.json")), error = function(e) NULL)
    if (!is.null(prior) && identical(prior$status, "complete")) {
      stop("Refusing to overwrite a completed MISTy output directory: ", output_dir)
    }
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_status("running")
  run_memory_gate(paste0("misty_", tolower(sample_id), "_read"))

  measurements_raw <- read_tsv(measurement_path, show_col_types = FALSE)
  positions <- read_tsv(position_path, show_col_types = FALSE)
  feature_manifest <- read_tsv(feature_manifest_path, show_col_types = FALSE)
  if (!identical(measurements_raw$barcode, positions$barcode)) {
    stop("Measurement and position barcodes differ")
  }
  measurements <- as.data.frame(measurements_raw[, -1, drop = FALSE])
  rownames(measurements) <- measurements_raw$barcode
  position_matrix <- as.data.frame(positions[, c("x_um", "y_um")])
  rownames(position_matrix) <- positions$barcode
  if (ncol(measurements) != 20L || nrow(measurements) != nrow(position_matrix)) {
    stop("Unexpected MISTy input dimensions")
  }
  if (!all(vapply(measurements, function(x) all(is.finite(x)), logical(1)))) {
    stop("Non-finite MISTy input")
  }
  if (max(abs(vapply(measurements, mean, numeric(1)))) > 1e-6 ||
      max(abs(vapply(measurements, sd, numeric(1)) - 1)) > 1e-3) {
    stop("MISTy input is not standardized as expected")
  }

  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  future::plan(future::multisession, workers = 4)
  on.exit(future::plan(future::sequential), add = TRUE)
  set.seed(seed)
  geometry_fallback_used <- FALSE
  geometry_jitter_epsilon_um <- 0
  position_matrix_model <- position_matrix
  initial_view <- create_initial_view(measurements, unique.id = paste0(sample_id, "_misty_v1"))
  views <- tryCatch(
    add_juxtaview(
      initial_view,
      position_matrix_model,
      neighbor.thr = 110,
      cached = FALSE,
      verbose = TRUE
    ),
    error = function(e) {
      if (!grepl("trifnd|deldir|clockwise", conditionMessage(e), ignore.case = TRUE)) {
        stop(e)
      }
      geometry_fallback_used <<- TRUE
      geometry_jitter_epsilon_um <<- 1e-3
      indices <- seq_len(nrow(position_matrix_model))
      position_matrix_model$x_um <<- position_matrix$x_um +
        geometry_jitter_epsilon_um * sin(indices * sqrt(2))
      position_matrix_model$y_um <<- position_matrix$y_um +
        geometry_jitter_epsilon_um * cos(indices * sqrt(3))
      message(
        "DELDIR_DEGENERACY_FALLBACK sample=", sample_id,
        " deterministic_coordinate_jitter_epsilon_um=", geometry_jitter_epsilon_um
      )
      retry_view <- create_initial_view(measurements, unique.id = paste0(sample_id, "_misty_v1"))
      add_juxtaview(
        retry_view,
        position_matrix_model,
        neighbor.thr = 110,
        cached = FALSE,
        verbose = TRUE
      )
    }
  )
  views <- add_paraview(
    views,
    position_matrix_model,
    l = 200,
    zoi = 110,
    family = "gaussian",
    cached = FALSE,
    verbose = TRUE
  )
  run_memory_gate(paste0("misty_", tolower(sample_id), "_model"))
  run_misty(
    views,
    results.folder = raw_dir,
    seed = seed,
    cv.folds = 10,
    cached = FALSE,
    num.trees = 100,
    num.threads = 1
  )
  result <- collect_results(raw_dir)

  improvements <- result$improvements
  contributions <- result$contributions
  importances <- result$importances
  improvement_stats <- result$improvements.stats
  contribution_stats <- result$contributions.stats
  importance_aggregated <- result$importances.aggregated
  improvements$sample <- sample_id
  contributions$sample <- sample_id
  importances$sample <- sample_id
  importance_aggregated$sample <- sample_id
  write_tsv(improvements, file.path(output_dir, paste0(sample_id, "_misty_improvements_v1.tsv")))
  write_tsv(contributions, file.path(output_dir, paste0(sample_id, "_misty_contributions_v1.tsv")))
  write_tsv(importances, file.path(output_dir, paste0(sample_id, "_misty_importances_v1.tsv")))
  write_tsv(improvement_stats, file.path(output_dir, paste0(sample_id, "_misty_improvement_stats_v1.tsv")))
  write_tsv(contribution_stats, file.path(output_dir, paste0(sample_id, "_misty_contribution_stats_v1.tsv")))
  write_tsv(importance_aggregated, file.path(output_dir, paste0(sample_id, "_misty_importances_aggregated_v1.tsv")))

  gains <- improvement_stats %>%
    filter(measure %in% c("gain.R2", "gain.RMSE")) %>%
    mutate(
      target_label = factor(safe_label(target, feature_manifest),
        levels = rev(safe_label(unique(target), feature_manifest))),
      measure = recode(measure, `gain.R2` = "R2 gain (percentage points)",
        `gain.RMSE` = "RMSE improvement (%)")
    )
  p_gain <- ggplot(gains, aes(x = mean, y = target_label, color = measure)) +
    geom_vline(xintercept = 0, color = "grey70", linewidth = 0.35) +
    geom_point(size = 2.1) +
    facet_wrap(~measure, scales = "free_x", ncol = 2) +
    labs(title = paste0(sample_id, ": MISTy multiview performance gain"), x = NULL, y = NULL) +
    theme_bw(base_size = 9) +
    theme(legend.position = "none", strip.background = element_rect(fill = "grey95"))
  ggsave(file.path(output_dir, paste0(sample_id, "_misty_performance_gains_v1.png")), p_gain,
    width = 11, height = 7.5, dpi = 300)
  ggsave(file.path(output_dir, paste0(sample_id, "_misty_performance_gains_v1.pdf")), p_gain,
    width = 11, height = 7.5)

  view_data <- contribution_stats %>%
    filter(view %in% c("intra", "juxta.110", "para.200")) %>%
    mutate(
      target_label = factor(safe_label(target, feature_manifest),
        levels = rev(safe_label(unique(target), feature_manifest))),
      view = factor(view, levels = c("intra", "juxta.110", "para.200"),
        labels = c("Intra", "Juxta <=110 um", "Para >110 um (l=200 um)"))
    )
  p_contrib <- ggplot(view_data, aes(x = fraction, y = target_label, fill = view)) +
    geom_col(width = 0.72) +
    scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
    scale_fill_manual(values = c("#4C78A8", "#F58518", "#54A24B")) +
    labs(title = paste0(sample_id, ": relative MISTy view contributions"),
      x = "Absolute meta-model coefficient fraction", y = NULL, fill = "View") +
    theme_bw(base_size = 9) +
    theme(legend.position = "top")
  ggsave(file.path(output_dir, paste0(sample_id, "_misty_view_contributions_v1.png")), p_contrib,
    width = 10, height = 7.5, dpi = 300)
  ggsave(file.path(output_dir, paste0(sample_id, "_misty_view_contributions_v1.pdf")), p_contrib,
    width = 10, height = 7.5)

  heat <- importances %>%
    filter(view %in% c("juxta.110", "para.200"), Predictor != Target) %>%
    mutate(
      predictor_label = safe_label(Predictor, feature_manifest),
      target_label = safe_label(Target, feature_manifest),
      view = factor(view, levels = c("juxta.110", "para.200"),
        labels = c("Juxta <=110 um", "Para >110 um (l=200 um)"))
    )
  cap <- max(abs(stats::quantile(heat$Importance, c(0.02, 0.98), na.rm = TRUE)))
  p_heat <- ggplot(heat, aes(x = predictor_label, y = target_label, fill = Importance)) +
    geom_tile() +
    facet_wrap(~view, ncol = 2) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426",
      midpoint = 0, limits = c(-cap, cap), oob = scales::squish) +
    labs(title = paste0(sample_id, ": MISTy relative predictor importance"),
      subtitle = "Scaled within each target-view and weighted by view support; not an effect direction",
      x = "Predictor", y = "Target", fill = "Importance") +
    theme_bw(base_size = 7.5) +
    theme(axis.text.x = element_text(angle = 55, hjust = 1), strip.background = element_rect(fill = "grey95"))
  ggsave(file.path(output_dir, paste0(sample_id, "_misty_juxta_para_importance_v1.png")), p_heat,
    width = 17, height = 9.5, dpi = 300)
  ggsave(file.path(output_dir, paste0(sample_id, "_misty_juxta_para_importance_v1.pdf")), p_heat,
    width = 17, height = 9.5)

  run_metadata <- list(
    status = "complete",
    sample = sample_id,
    n_spots = nrow(measurements),
    n_features = ncol(measurements),
    feature_ids = colnames(measurements),
    views = list(
      intra = "same spot",
      juxta = list(neighbor_threshold_um = 110, geometry = "Delaunay neighbors within threshold"),
      para = list(family = "gaussian", l_um = 200, zoi_um = 110)
    ),
    random_forest = list(num_trees = 100, num_threads_per_model = 1),
    future_workers = 4,
    cv_folds = 10,
    seed = seed,
    mistyR_version = as.character(packageVersion("mistyR")),
    interpretation = paste(
      "MISTy importances are multivariate predictive associations; scaled importance is not an effect",
      "direction and must not be interpreted as a linear or causal relationship."
    ),
    measurement_path = measurement_path,
    measurement_md5 = md5sum_file(measurement_path),
    position_path = position_path,
    position_md5 = md5sum_file(position_path),
    geometry_fallback = list(
      used = geometry_fallback_used,
      reason = if (geometry_fallback_used) "deldir degenerate regular-lattice triangulation" else "not needed",
      deterministic_coordinate_jitter_epsilon_um = geometry_jitter_epsilon_um,
      maximum_possible_axis_shift_relative_to_110um = geometry_jitter_epsilon_um / 110
    ),
    completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
  )
  write_json(run_metadata, file.path(output_dir, "run_metadata_v1.json"), pretty = TRUE, auto_unbox = TRUE)
  write_status("complete")
  message("MISTY_SAMPLE_ALL_PASS sample=", sample_id, " spots=", nrow(measurements), " features=", ncol(measurements))
}

tryCatch(
  main(),
  error = function(e) {
    write_status("failed", conditionMessage(e))
    stop(e)
  }
)
