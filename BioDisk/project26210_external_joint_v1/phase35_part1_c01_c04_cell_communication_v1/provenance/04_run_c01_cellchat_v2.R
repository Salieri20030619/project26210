suppressPackageStartupMessages({
  library(CellChat)
  library(Matrix)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) stop("Usage: Rscript script input_dir db.rds output_dir [comma_separated_patients]")
input_dir <- args[[1]]
db_path <- args[[2]]
out <- args[[3]]
patient_selector <- if (length(args) == 4L) args[[4]] else ""
if (dir.exists(out)) stop("Refusing to overwrite output directory")
dir.create(out, recursive = TRUE)
dir.create(file.path(out, "patient_rds"))

expr <- readMM(gzfile(file.path(input_dir, "cellchat_expression_genes_by_cells_v1.mtx.gz")))
genes <- read.delim(file.path(input_dir, "cellchat_genes_v1.tsv"), check.names = FALSE, stringsAsFactors = FALSE)$gene
meta <- read.delim(file.path(input_dir, "cellchat_cell_metadata_v1.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
if (nrow(expr) != length(genes) || ncol(expr) != nrow(meta)) stop("Expression/metadata dimension mismatch")
rownames(expr) <- genes
colnames(expr) <- meta$exact_cell_id
rownames(meta) <- meta$exact_cell_id
db <- readRDS(db_path)
future::plan("sequential")

group_counts <- as.data.frame(table(meta$patient, meta$communication_group), stringsAsFactors = FALSE)
names(group_counts) <- c("patient", "communication_group", "n_cells")
eligible <- unique(meta$patient)
eligible <- eligible[vapply(eligible, function(patient) {
  z <- group_counts[group_counts$patient == patient, ]
  high <- z$n_cells[z$communication_group == "TYMP-high Malignant"]
  low <- z$n_cells[z$communication_group == "TYMP-low Malignant"]
  length(high) == 1L && length(low) == 1L && high >= 20L && low >= 20L
}, logical(1))]
eligible <- sort(eligible)
if (nzchar(patient_selector)) {
  requested <- strsplit(patient_selector, ",", fixed = TRUE)[[1]]
  missing_requested <- setdiff(requested, eligible)
  if (length(missing_requested)) stop("Requested patients are not eligible: ", paste(missing_requested, collapse = ","))
  eligible <- requested
}
if (!length(eligible)) stop("No eligible patients")

all_interactions <- list()
all_endpoints <- list()
all_lr_diffs <- list()
run_qc <- list()
macrophages <- c("SPP1-state Macrophage", "C1QC-state Macrophage", "Other Macrophage")

edge_sum <- function(tab, sources, targets) {
  value <- tab$prob[tab$source %in% sources & tab$target %in% targets]
  if (!length(value)) return(0)
  sum(value, na.rm = TRUE)
}

for (patient in eligible) {
  message(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), " ", patient)
  idx <- which(meta$patient == patient)
  pmeta <- meta[idx, , drop = FALSE]
  present <- unique(pmeta$communication_group)
  cellchat <- createCellChat(object = expr[, idx, drop = FALSE], meta = pmeta, group.by = "communication_group")
  cellchat@DB <- db
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                                population.size = FALSE, raw.use = TRUE)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)
  tab <- subsetCommunication(cellchat, thresh = 1)
  if (is.null(tab) || !nrow(tab)) {
    tab <- data.frame(source = character(), target = character(), ligand = character(),
                      receptor = character(), interaction_name = character(), pathway_name = character(),
                      prob = numeric(), pval = numeric())
  }
  tab$patient <- patient
  tab$study <- unique(pmeta$study_v1)[[1]]
  all_interactions[[patient]] <- tab

  endpoint_specs <- list(
    malignant_to_macrophage = list(high = "TYMP-high Malignant", low = "TYMP-low Malignant", target = macrophages),
    macrophage_to_malignant = list(high = macrophages, low = macrophages, target = c("TYMP-high Malignant", "TYMP-low Malignant")),
    malignant_to_endothelial = list(high = "TYMP-high Malignant", low = "TYMP-low Malignant", target = "Endothelial"),
    endothelial_to_malignant = list(high = "Endothelial", low = "Endothelial", target = c("TYMP-high Malignant", "TYMP-low Malignant"))
  )
  for (endpoint in names(endpoint_specs)) {
    spec <- endpoint_specs[[endpoint]]
    if (endpoint %in% c("malignant_to_macrophage", "malignant_to_endothelial")) {
      target_present <- intersect(spec$target, present)
      eligible_endpoint <- length(target_present) > 0L
      high_score <- if (eligible_endpoint) edge_sum(tab, spec$high, target_present) else NA_real_
      low_score <- if (eligible_endpoint) edge_sum(tab, spec$low, target_present) else NA_real_
    } else {
      source_present <- intersect(spec$high, present)
      eligible_endpoint <- length(source_present) > 0L
      high_score <- if (eligible_endpoint) edge_sum(tab, source_present, "TYMP-high Malignant") else NA_real_
      low_score <- if (eligible_endpoint) edge_sum(tab, source_present, "TYMP-low Malignant") else NA_real_
    }
    all_endpoints[[paste(patient, endpoint)]] <- data.frame(
      study = unique(pmeta$study_v1)[[1]], patient = patient, endpoint = endpoint,
      eligible = eligible_endpoint, high_score = high_score, low_score = low_score,
      high_minus_low = high_score - low_score
    )
  }

  if (nrow(tab)) {
    pair_id <- if ("interaction_name" %in% names(tab)) tab$interaction_name else paste(tab$ligand, tab$receptor, sep = "_")
    tab$pair_id <- pair_id
    for (compartment in c("macrophage", "endothelial")) {
      targets <- if (compartment == "macrophage") intersect(macrophages, present) else intersect("Endothelial", present)
      if (!length(targets)) next
      outgoing <- aggregate(prob ~ pair_id + pathway_name, tab[tab$source %in% c("TYMP-high Malignant", "TYMP-low Malignant") & tab$target %in% targets, ], sum)
      if (nrow(outgoing)) {
        high_tab <- aggregate(prob ~ pair_id + pathway_name, tab[tab$source == "TYMP-high Malignant" & tab$target %in% targets, ], sum)
        low_tab <- aggregate(prob ~ pair_id + pathway_name, tab[tab$source == "TYMP-low Malignant" & tab$target %in% targets, ], sum)
        merged <- merge(high_tab, low_tab, by = c("pair_id", "pathway_name"), all = TRUE, suffixes = c("_high", "_low"))
        merged[is.na(merged)] <- 0
        merged$direction_family <- paste0("malignant_to_", compartment)
        merged$patient <- patient
        merged$study <- unique(pmeta$study_v1)[[1]]
        merged$high_minus_low <- merged$prob_high - merged$prob_low
        all_lr_diffs[[paste(patient, compartment, "out")]] <- merged
      }
    }
  }
  saveRDS(cellchat, file.path(out, "patient_rds", paste0(gsub("[^A-Za-z0-9]+", "_", patient), "_cellchat_v1.rds")))
  run_qc[[patient]] <- data.frame(study = unique(pmeta$study_v1)[[1]], patient = patient,
                                  n_cells = length(idx), n_groups = length(present),
                                  n_inferred_rows = nrow(tab), status = "PASS")
  rm(cellchat); gc()
}

interactions <- do.call(rbind, all_interactions)
endpoints <- do.call(rbind, all_endpoints)
lr_diffs <- if (length(all_lr_diffs)) do.call(rbind, all_lr_diffs) else data.frame()
qc <- do.call(rbind, run_qc)

summaries <- lapply(split(endpoints[endpoints$eligible & is.finite(endpoints$high_minus_low), ], endpoints$endpoint), function(z) {
  p <- if (nrow(z) >= 3L && any(z$high_minus_low != 0)) wilcox.test(z$high_minus_low, mu = 0, alternative = "two.sided", exact = FALSE)$p.value else 1
  data.frame(endpoint = z$endpoint[[1]], n_eligible_patients = nrow(z),
             n_positive = sum(z$high_minus_low > 0), n_negative = sum(z$high_minus_low < 0), n_zero = sum(z$high_minus_low == 0),
             median_high_minus_low = median(z$high_minus_low), wilcoxon_p_two_sided = p)
})
summary <- do.call(rbind, summaries)
summary$BH_FDR_4_endpoints <- p.adjust(summary$wilcoxon_p_two_sided, method = "BH")
summary$statistical_support <- ifelse(summary$n_eligible_patients >= 10, summary$BH_FDR_4_endpoints < 0.05,
                                      ifelse(summary$n_eligible_patients >= 3, summary$wilcoxon_p_two_sided < 0.05, FALSE))
summary$direction <- ifelse(summary$median_high_minus_low > 0, "positive", ifelse(summary$median_high_minus_low < 0, "negative", "zero"))

write_gz <- function(x, path) {
  part <- paste0(path, ".part")
  con <- gzfile(part, "wt")
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE)
  close(con)
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}
write_plain <- function(x, path) {
  part <- paste0(path, ".part")
  write.table(x, part, sep = "\t", quote = FALSE, row.names = FALSE)
  if (!file.rename(part, path)) stop("Atomic rename failed: ", path)
}
write_gz(interactions, file.path(out, "c01_cellchat_patient_interactions_v1.tsv.gz"))
write_plain(endpoints, file.path(out, "c01_cellchat_patient_endpoint_effects_v1.tsv"))
write_plain(summary, file.path(out, "c01_cellchat_endpoint_summary_v1.tsv"))
write_gz(lr_diffs, file.path(out, "c01_cellchat_lr_pair_differences_v1.tsv.gz"))
write_plain(qc, file.path(out, "c01_cellchat_patient_run_qc_v1.tsv"))
payload <- list(status = "PASS_C01_CELLCHAT_ANALYSIS_COMPLETE_PRESENTATION_PENDING",
                completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
                n_patients = length(eligible), n_interaction_rows = nrow(interactions),
                two_sided_testing = TRUE, direction_prespecified = FALSE, figures_created = FALSE)
part <- file.path(out, "status.json.part")
write_json(payload, part, pretty = TRUE, auto_unbox = TRUE)
if (!file.rename(part, file.path(out, "status.json"))) stop("Atomic status rename failed")
