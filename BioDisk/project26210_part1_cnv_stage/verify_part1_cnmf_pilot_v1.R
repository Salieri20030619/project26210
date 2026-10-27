#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(digest)
})

PROJECT <- "/mnt/BioDisk/ludingyu/project26210"
ROOT <- file.path(PROJECT, "results/part1_rebuild_v1")
INPUT <- file.path(ROOT, "05_cnmf/input_v1")
PILOT <- file.path(ROOT, "05_cnmf/pilot_v1")
OUT <- file.path(ROOT, "05_cnmf/pilot_verification_v1")
MEMORY_GATE <- file.path(PROJECT, "check_memory_gate.sh")
RANKS <- c(4L, 6L, 8L)
N_RUNS <- 3L
EXPECTED_MATRIX_SHA256 <- "65aeb7e54d6f5c1f526882b71ddcb0f96801dde392d4d58f593c6a50e59f2935"

if (dir.exists(OUT) && length(list.files(OUT, all.files = TRUE, no.. = TRUE)) > 0L) {
  stop("Refusing to overwrite non-empty pilot verification directory: ", OUT)
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
checks <- list()
check <- function(name, passed, evidence) {
  checks[[length(checks) + 1L]] <<- data.table(
    check = name,
    status = if (isTRUE(passed)) "PASS" else "FAIL",
    evidence = paste(evidence, collapse = ";")
  )
  current <- rbindlist(checks, fill = TRUE)
  fwrite(current, file.path(OUT, "verification_checks_v1.tsv"), sep = "\t")
  if (!isTRUE(passed)) stop("Verification failed: ", name, ": ", paste(evidence, collapse = ";"))
}

gate_output <- system2(MEMORY_GATE, "cnmf_pilot_independent_verification_v1", stdout = TRUE, stderr = TRUE)
gate_status <- attr(gate_output, "status")
if (is.null(gate_status)) gate_status <- 0L
gate_text <- paste(gate_output, collapse = " ")
message(gate_text)
check("memory_gate", gate_status == 0L && grepl("PASS", gate_text, fixed = TRUE), gate_text)
check("pilot_complete", file.exists(file.path(PILOT, "RUN_COMPLETE")), file.path(PILOT, "RUN_COMPLETE"))
check("input_verification_complete", file.exists(file.path(ROOT, "05_cnmf/input_verification_v1/RUN_COMPLETE")), "terminal")

input_manifest <- fread(file.path(INPUT, "file_manifest_v1.tsv"), colClasses = "character")
matrix_row <- input_manifest[file == "cnmf_hvg_scaled_genes_by_cells_v1.mtx.gz"]
matrix_path <- file.path(INPUT, matrix_row$file)
matrix_sha <- digest(matrix_path, algo = "sha256", file = TRUE)
check("input_matrix_manifest_unique", nrow(matrix_row) == 1L, nrow(matrix_row))
check("input_matrix_hash", identical(matrix_sha, EXPECTED_MATRIX_SHA256) && identical(matrix_sha, matrix_row$sha256), matrix_sha)

metadata <- fromJSON(file.path(PILOT, "run_metadata_v1.json"), simplifyVector = TRUE)
check("metadata_schema", identical(metadata$schema_version, "part1_cnmf_pilot_v1"), metadata$schema_version)
check("metadata_purpose", grepl("not a biological rank selection", metadata$purpose, fixed = TRUE), metadata$purpose)
check("metadata_input_hash", identical(metadata$input_matrix_sha256, matrix_sha), metadata$input_matrix_sha256)
check("metadata_dimensions", metadata$n_genes == 2500L && metadata$n_cells == 11950L, c(metadata$n_genes, metadata$n_cells))
check("metadata_ranks", identical(as.integer(metadata$ranks), RANKS), metadata$ranks)
check("metadata_parameters", metadata$n_runs_per_rank == N_RUNS && metadata$maxit == 75L && metadata$seed_base == 26210L, c(metadata$n_runs_per_rank, metadata$maxit, metadata$seed_base))
check("metadata_environment", metadata$rcppml_threads == 4L && identical(metadata$RcppML_version, "0.3.7.1"), c(metadata$rcppml_threads, metadata$RcppML_version))
check("metadata_no_dependency_changes", identical(metadata$dependency_changes, FALSE), metadata$dependency_changes)

gates <- fread(file.path(PILOT, "memory_gate_log_v1.tsv"))
expected_stages <- c("pilot_load_input_v1",
                     paste0("pilot_rank4_run", 1:3), "pilot_rank4_consensus",
                     paste0("pilot_rank6_run", 1:3), "pilot_rank6_consensus",
                     paste0("pilot_rank8_run", 1:3), "pilot_rank8_consensus")
check("memory_gate_count", nrow(gates) == 13L, nrow(gates))
check("memory_gate_stages", identical(gates$stage, expected_stages), gates$stage)
check("memory_gates_all_pass", all(gates$exit_status == 0L) && all(grepl("PASS", gates$gate_output, fixed = TRUE)), sum(grepl("PASS", gates$gate_output, fixed = TRUE)))

run_qc <- fread(file.path(PILOT, "run_qc_v1.tsv"))
expected_runs <- CJ(rank = RANKS, run = seq_len(N_RUNS), sorted = TRUE)
expected_runs[, seed := as.integer(26210L + rank * 100L + run)]
setorder(run_qc, rank, run)
check("run_qc_rows", nrow(run_qc) == 9L && !anyDuplicated(run_qc[, .(rank, run)]), nrow(run_qc))
check("run_qc_grid", all(run_qc$rank == expected_runs$rank) && all(run_qc$run == expected_runs$run) && all(run_qc$seed == expected_runs$seed), "ranks 4/6/8 x runs 1/2/3")
check("run_qc_iterations", all(run_qc$iter >= 1 & run_qc$iter <= 75), range(run_qc$iter))
check("run_qc_tolerance", all(is.finite(run_qc$tolerance)) && all(run_qc$tolerance >= 0), range(run_qc$tolerance))
check("run_qc_mse", all(is.finite(run_qc$mse)) && all(run_qc$mse > 0), range(run_qc$mse))
check("run_qc_elapsed", all(is.finite(run_qc$elapsed_seconds)) && all(run_qc$elapsed_seconds > 0), range(run_qc$elapsed_seconds))
check("run_qc_fit_valid", all(run_qc$fit_finite_nonnegative), sum(run_qc$fit_finite_nonnegative))

component <- fread(file.path(PILOT, "component_assignments_v1.tsv"))
check("component_rows", nrow(component) == sum(RANKS) * N_RUNS && !anyDuplicated(component[, .(rank, component_id)]), nrow(component))
check("component_rank_counts", identical(component[, .N, by = rank][order(rank)]$N, RANKS * N_RUNS), component[, .N, by = rank]$N)
check("component_silhouette_finite", all(is.finite(component$silhouette)) && all(component$silhouette >= -1) && all(component$silhouette <= 1), range(component$silhouette))
component_cluster <- component[, .(n_components = .N, n_unique_runs = uniqueN(run)), by = .(rank, consensus_cluster)]
check("component_cluster_counts", identical(component_cluster[, .N, by = rank][order(rank)]$N, RANKS), component_cluster[, .N, by = rank]$N)
check("component_all_runs_per_cluster", all(component_cluster$n_unique_runs == N_RUNS), range(component_cluster$n_unique_runs))
check("component_three_per_cluster", all(component_cluster$n_components == N_RUNS), range(component_cluster$n_components))

rank_qc <- fread(file.path(PILOT, "rank_qc_v1.tsv"))
setorder(rank_qc, rank)
check("rank_qc_rows", identical(rank_qc$rank, RANKS), rank_qc$rank)
check("rank_qc_runs", all(rank_qc$n_runs == N_RUNS), rank_qc$n_runs)
check("rank_qc_full_run_fraction", all(abs(rank_qc$min_run_fraction - 1) < 1e-12), rank_qc$min_run_fraction)
check("rank_qc_mean_silhouette_gate", all(rank_qc$mean_component_silhouette >= 0.70), rank_qc$mean_component_silhouette)
recomputed_sil <- component[, .(value = mean(silhouette)), by = rank][order(rank)]$value
check("rank_qc_silhouette_reconciles", max(abs(rank_qc$mean_component_silhouette - recomputed_sil)) < 1e-12, max(abs(rank_qc$mean_component_silhouette - recomputed_sil)))
recomputed_mse <- run_qc[, .(value = median(mse)), by = rank][order(rank)]$value
check("rank_qc_mse_reconciles", max(abs(rank_qc$median_run_mse - recomputed_mse)) < 1e-12, max(abs(rank_qc$median_run_mse - recomputed_mse)))

cells <- fread(file.path(INPUT, "cnmf_cells_v1.tsv.gz"), colClasses = c(exact_cell_id = "character"))$exact_cell_id
genes <- fread(file.path(INPUT, "cnmf_hvg_genes_v1.tsv"), colClasses = c(gene = "character"))$gene
models <- readRDS(file.path(PILOT, "pilot_models_v1.rds"))
check("model_parameters", identical(as.integer(models$ranks), RANKS) && models$n_runs == N_RUNS && models$maxit == 75L && identical(models$input_matrix_sha256, matrix_sha), "exact")
check("model_rank_names", identical(names(models$fits), as.character(RANKS)) && identical(names(models$consensus), as.character(RANKS)), names(models$fits))

for (rank_value in RANKS) {
  rank_name <- as.character(rank_value)
  fits <- models$fits[[rank_name]]
  result <- models$consensus[[rank_name]]
  check(paste0("rank", rank_value, "_fit_count"), length(fits) == N_RUNS, length(fits))
  fit_valid <- all(vapply(fits, function(fit) {
    identical(dim(fit$w), c(2500L, rank_value)) &&
      identical(dim(fit$h), c(rank_value, 11950L)) && length(fit$d) == rank_value &&
      all(is.finite(fit$w)) && all(fit$w >= 0) && all(is.finite(fit$h)) &&
      all(fit$h >= 0) && all(is.finite(fit$d)) && all(fit$d >= 0)
  }, logical(1)))
  check(paste0("rank", rank_value, "_fits_valid"), fit_valid, fit_valid)

  usage <- fread(file.path(PILOT, paste0("consensus_usage_rank", rank_value, "_v1.tsv.gz")))
  loading <- fread(file.path(PILOT, paste0("consensus_loading_rank", rank_value, "_v1.tsv.gz")))
  usage_values <- as.matrix(usage[, setdiff(names(usage), "exact_cell_id"), with = FALSE])
  loading_values <- as.matrix(loading[, setdiff(names(loading), "gene"), with = FALSE])
  check(paste0("rank", rank_value, "_usage_dimensions"), identical(dim(usage_values), c(11950L, rank_value)), dim(usage_values))
  check(paste0("rank", rank_value, "_usage_cell_order"), identical(usage$exact_cell_id, cells), "exact")
  check(paste0("rank", rank_value, "_usage_values"), all(is.finite(usage_values)) && all(usage_values >= 0) && max(abs(rowSums(usage_values) - 1)) < 1e-12, max(abs(rowSums(usage_values) - 1)))
  check(paste0("rank", rank_value, "_loading_dimensions"), identical(dim(loading_values), c(2500L, rank_value)), dim(loading_values))
  check(paste0("rank", rank_value, "_loading_gene_order"), identical(loading$gene, genes), "exact")
  check(paste0("rank", rank_value, "_loading_values"), all(is.finite(loading_values)) && all(loading_values >= 0) && max(abs(colSums(loading_values) - 1)) < 1e-12, max(abs(colSums(loading_values) - 1)))
  check(paste0("rank", rank_value, "_usage_matches_rds"), max(abs(usage_values - result$usage)) < 1e-12, max(abs(usage_values - result$usage)))
  check(paste0("rank", rank_value, "_loading_matches_rds"), max(abs(loading_values - result$loading)) < 1e-12, max(abs(loading_values - result$loading)))
}

manifest <- fread(file.path(PILOT, "file_manifest_v1.tsv"), colClasses = "character")
check("manifest_unique", !anyDuplicated(manifest$file), uniqueN(manifest$file))
for (i in seq_len(nrow(manifest))) {
  path <- file.path(PILOT, manifest$file[i])
  valid <- file.exists(path) && file.info(path)$size == as.numeric(manifest$size_bytes[i]) &&
    identical(digest(path, algo = "sha256", file = TRUE), manifest$sha256[i])
  check(paste0("hash_", manifest$file[i]), valid, manifest$sha256[i])
}

final_checks <- rbindlist(checks, fill = TRUE)
run_metadata <- list(
  schema_version = "part1_cnmf_pilot_verification_v1",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  n_checks = nrow(final_checks),
  n_pass = sum(final_checks$status == "PASS"),
  n_fail = sum(final_checks$status == "FAIL"),
  computational_gate = "PASS",
  biological_rank_selected = FALSE,
  dependency_changes = FALSE
)
write_json(run_metadata, file.path(OUT, "run_metadata_v1.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(OUT, "RUN_COMPLETE"))
print(run_metadata)
cat("PART1_CNMF_PILOT_VERIFICATION_V1_PASS\n")
