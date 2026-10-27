#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

set.seed(26210)
project <- "/mnt/BioDisk/ludingyu/project26210"
root <- file.path(project, "results", "part1_rebuild_v1")
input_dir <- file.path(root, "07_pyscenic_v1", "input")
cnmf_dir <- file.path(root, "05_cnmf")
clone_dir <- file.path(root, "04_subclones_tymp", "cnv_clone_tymp_cnmf_association_v1")
out <- file.path(cnmf_dir, "cnv_clone_program_function_v1")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)

gate_status <- system2(
  file.path(project, "check_memory_gate.sh"),
  "part1_cnv_cnmf_function_v1",
  stdout = file.path(out, "memory_gate_v1.txt"),
  stderr = file.path(out, "memory_gate_v1.txt")
)
if (!identical(gate_status, 0L)) stop("Memory gate refused CNV/cNMF function analysis")

required <- c(
  file.path(input_dir, "RUN_COMPLETE"),
  file.path(input_dir, "raw_counts_genes_by_cells_v1.mtx.gz"),
  file.path(input_dir, "genes_v1.tsv"),
  file.path(input_dir, "cells_v1.tsv"),
  file.path(input_dir, "cell_metadata_v1.tsv.gz"),
  file.path(cnmf_dir, "full_ranks4_8_v1", "consensus_loading_rank7_v1.tsv.gz"),
  file.path(clone_dir, "patient_factor_high_clone_tymp_contrasts_v1.tsv")
)
stopifnot(all(file.exists(required)))

genes <- fread(file.path(input_dir, "genes_v1.tsv"))$gene
cells <- fread(file.path(input_dir, "cells_v1.tsv"))$exact_cell_id
meta <- fread(file.path(input_dir, "cell_metadata_v1.tsv.gz"))
stopifnot(identical(meta$exact_cell_id, cells), !anyDuplicated(cells), !anyDuplicated(genes))

con <- gzfile(file.path(input_dir, "raw_counts_genes_by_cells_v1.mtx.gz"), "rb")
counts <- readMM(con)
close(con)
counts <- as(counts, "dgCMatrix")
stopifnot(identical(dim(counts), c(length(genes), length(cells))))
rownames(counts) <- genes
colnames(counts) <- cells
stopifnot(length(counts@x) == 35713877L, all(is.finite(counts@x)), all(counts@x >= 0), all(counts@x == round(counts@x)))

contrasts <- fread(file.path(clone_dir, "patient_factor_high_clone_tymp_contrasts_v1.tsv"))
contrasts <- contrasts[distance_threshold == 0.25 & factor == "F1"]
stopifnot(nrow(contrasts) == 14L, !anyDuplicated(contrasts$infercnv_unit))
clone_map <- contrasts[, .(curated_patient_id = infercnv_unit, factor_high_clone_id)]
meta <- merge(meta, clone_map, by = "curated_patient_id", all.x = TRUE, sort = FALSE)
meta[, order_index__ := match(exact_cell_id, cells)]
setorder(meta, order_index__)
meta[, order_index__ := NULL]
stopifnot(identical(meta$exact_cell_id, cells))
meta[, comparison_group := fifelse(
  clone_size_class_v1 == "primary_clone_ge20_v1" & formal_inferred_cnv_clone_v1 == factor_high_clone_id,
  "F1_high_primary_clone",
  fifelse(clone_size_class_v1 == "primary_clone_ge20_v1" & !is.na(factor_high_clone_id), "other_primary_clones", NA_character_)
)]

selected <- which(!is.na(meta$comparison_group))
sample_meta <- unique(meta[selected, .(study, curated_patient_id, comparison_group)])
sample_meta[, sample_key := paste(curated_patient_id, comparison_group, sep = "||")]
sample_meta[, group_order__ := match(comparison_group, c("other_primary_clones", "F1_high_primary_clone"))]
setorder(sample_meta, curated_patient_id, group_order__)
sample_meta[, group_order__ := NULL]
stopifnot(nrow(sample_meta) == 28L, all(table(sample_meta$curated_patient_id) == 2L))
cell_key <- paste(meta$curated_patient_id[selected], meta$comparison_group[selected], sep = "||")
aggregate_map <- sparseMatrix(
  i = seq_along(selected),
  j = match(cell_key, sample_meta$sample_key),
  x = 1,
  dims = c(length(selected), nrow(sample_meta))
)
pseudobulk <- counts[, selected, drop = FALSE] %*% aggregate_map
colnames(pseudobulk) <- sample_meta$sample_key
rownames(pseudobulk) <- genes
sample_cell_counts <- tabulate(match(cell_key, sample_meta$sample_key), nbins = nrow(sample_meta))
sample_meta[, n_cells := sample_cell_counts]
sample_meta[, library_size := as.numeric(Matrix::colSums(pseudobulk))]
fwrite(sample_meta, file.path(out, "pseudobulk_sample_qc_v1.tsv"), sep = "\t")

group <- factor(sample_meta$comparison_group, levels = c("other_primary_clones", "F1_high_primary_clone"))
patient <- factor(sample_meta$curated_patient_id)
design <- model.matrix(~ patient + group)
y <- DGEList(counts = as.matrix(pseudobulk), genes = data.frame(gene = genes))
keep <- filterByExpr(y, design = design, min.count = 10)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
qlf <- glmQLFTest(fit, coef = "groupF1_high_primary_clone")
de_frame <- topTags(qlf, n = Inf, sort.by = "PValue")$table
de <- if ("gene" %in% names(de_frame)) as.data.table(de_frame) else as.data.table(de_frame, keep.rownames = "gene")
setnames(de, c("logFC", "logCPM", "F", "PValue", "FDR"), c("edgeR_logFC", "edgeR_logCPM", "edgeR_QLF", "edgeR_PValue", "edgeR_FDR"))

logcpm <- cpm(y, log = TRUE, prior.count = 2)
patient_levels <- levels(patient)
deltas <- sapply(patient_levels, function(p) {
  hi <- which(sample_meta$curated_patient_id == p & sample_meta$comparison_group == "F1_high_primary_clone")
  lo <- which(sample_meta$curated_patient_id == p & sample_meta$comparison_group == "other_primary_clones")
  logcpm[, hi] - logcpm[, lo]
})
colnames(deltas) <- patient_levels
paired_p <- apply(deltas, 1, function(x) {
  if (all(abs(x) < 1e-12)) return(1)
  suppressWarnings(wilcox.test(x, mu = 0, exact = FALSE)$p.value)
})
paired <- data.table(
  gene = rownames(deltas),
  paired_median_logCPM_delta = apply(deltas, 1, median),
  paired_mean_logCPM_delta = rowMeans(deltas),
  n_patients_positive = rowSums(deltas > 0),
  n_patients_negative = rowSums(deltas < 0),
  paired_wilcox_PValue = paired_p,
  paired_wilcox_FDR = p.adjust(paired_p, "BH")
)
de <- merge(de, paired, by = "gene", all.x = TRUE, sort = FALSE)
setorder(de, edgeR_PValue)
fwrite(de, file.path(out, "F1_high_clone_vs_other_paired_pseudobulk_DE_v1.tsv.gz"), sep = "\t")

top_delta_genes <- unique(c(head(de[edgeR_logFC > 0]$gene, 75), head(de[edgeR_logFC < 0]$gene, 75)))
delta_long <- as.data.table(deltas[top_delta_genes, , drop = FALSE], keep.rownames = "gene")
delta_long <- melt(delta_long, id.vars = "gene", variable.name = "curated_patient_id", value.name = "logCPM_delta")
delta_long <- merge(delta_long, unique(sample_meta[, .(curated_patient_id, study)]), by = "curated_patient_id", all.x = TRUE)
fwrite(delta_long, file.path(out, "top_DE_patient_logCPM_deltas_v1.tsv.gz"), sep = "\t")

read_gmt <- function(path, collection) {
  lines <- readLines(path, warn = FALSE)
  result <- lapply(lines, function(line) {
    fields <- strsplit(line, "\t", fixed = FALSE)[[1]]
    geneset <- unique(fields[-c(1, 2)])
    list(name = paste(collection, fields[1], sep = "::"), genes = geneset[nzchar(geneset)])
  })
  setNames(lapply(result, `[[`, "genes"), vapply(result, `[[`, character(1), "name"))
}

hallmark_gmt <- file.path(project, "envs", "part2_6_2_7_20260726_v1", "sources", "VISION", "vignettes", "data", "h.all.v5.2.symbols.gmt")
kegg_gmt <- file.path(project, "envs", "part2_6_2_7_20260726_v1", "sources", "scMetabolism", "data", "KEGG_metabolism_nc.gmt")
reactome_gmt <- file.path(project, "envs", "part2_6_2_7_20260726_v1", "sources", "scMetabolism", "data", "REACTOME_metabolism.gmt")
stopifnot(file.exists(hallmark_gmt), file.exists(kegg_gmt), file.exists(reactome_gmt))
pathways <- c(read_gmt(hallmark_gmt, "Hallmark_v5.2"), read_gmt(kegg_gmt, "KEGG_metabolism"), read_gmt(reactome_gmt, "Reactome_metabolism"))

ora_one <- function(selected_genes, universe, label, direction) {
  selected_genes <- intersect(unique(selected_genes), universe)
  rbindlist(lapply(names(pathways), function(pathway) {
    members <- intersect(pathways[[pathway]], universe)
    overlap <- intersect(selected_genes, members)
    data.table(
      analysis = label,
      direction = direction,
      pathway = pathway,
      selected_size = length(selected_genes),
      pathway_size = length(members),
      overlap_size = length(overlap),
      overlap_genes = paste(overlap, collapse = ","),
      PValue = phyper(length(overlap) - 1, length(members), length(universe) - length(members), length(selected_genes), lower.tail = FALSE)
    )
  }))[pathway_size >= 5 & pathway_size <= 500][, FDR := p.adjust(PValue, "BH")]
}

loadings <- fread(file.path(cnmf_dir, "full_ranks4_8_v1", "consensus_loading_rank7_v1.tsv.gz"))
stopifnot(nrow(loadings) == 2500L, identical(names(loadings), c("gene", paste0("F", 1:7))))
program_ora <- rbindlist(lapply(paste0("F", 1:7), function(factor_name) {
  top <- loadings[order(-get(factor_name))][1:100, gene]
  ora_one(top, loadings$gene, factor_name, "top100_loading")
}))
fwrite(program_ora[order(analysis, FDR, PValue)], file.path(out, "rank7_program_pathway_ORA_v1.tsv"), sep = "\t")

de_universe <- de$gene
de_ora <- rbindlist(list(
  ora_one(head(de[edgeR_logFC > 0][order(edgeR_PValue)]$gene, 250), de_universe, "F1_high_clone_vs_other", "up"),
  ora_one(head(de[edgeR_logFC < 0][order(edgeR_PValue)]$gene, 250), de_universe, "F1_high_clone_vs_other", "down")
))
fwrite(de_ora[order(direction, FDR, PValue)], file.path(out, "F1_high_clone_DE_pathway_ORA_v1.tsv"), sep = "\t")

program_top <- rbindlist(lapply(paste0("F", 1:7), function(f) {
  ord <- order(-loadings[[f]])[1:100]
  data.table(factor = f, gene = loadings$gene[ord], gene_rank = 1:100, loading = loadings[[f]][ord])
}))
fwrite(program_top, file.path(out, "rank7_program_top100_genes_v1.tsv"), sep = "\t")

plot_program <- program_ora[overlap_size >= 2][order(FDR, PValue), head(.SD, 5), by = analysis]
plot_paths <- unique(plot_program$pathway)
plot_grid <- CJ(analysis = paste0("F", 1:7), pathway = plot_paths)
plot_grid <- merge(plot_grid, program_ora[, .(analysis, pathway, FDR, overlap_size)], by = c("analysis", "pathway"), all.x = TRUE)
plot_grid[is.na(FDR), `:=`(FDR = 1, overlap_size = 0L)]
plot_grid[, pathway_label := sub("^[^:]+::", "", pathway)]
p1 <- ggplot(plot_grid, aes(analysis, factor(pathway_label, levels = rev(unique(plot_grid$pathway_label))))) +
  geom_point(aes(size = overlap_size, colour = pmin(-log10(pmax(FDR, 1e-12)), 12))) +
  scale_colour_viridis_c(name = "-log10 FDR") +
  scale_size_continuous(name = "Overlap") +
  labs(x = "rank-7 cNMF program", y = NULL, title = "cNMF program pathway annotation") +
  theme_bw(base_size = 8) + theme(panel.grid.minor = element_blank(), axis.text.y = element_text(size = 6))
ggsave(file.path(out, "figures", "rank7_program_pathway_ORA_v1.png"), p1, width = 8.5, height = max(5, length(plot_paths) * 0.18), dpi = 300, bg = "white")
ggsave(file.path(out, "figures", "rank7_program_pathway_ORA_v1.pdf"), p1, width = 8.5, height = max(5, length(plot_paths) * 0.18))

volcano <- copy(de)
volcano[, neglog10FDR := -log10(pmax(edgeR_FDR, 1e-300))]
volcano[, class := fifelse(edgeR_FDR < 0.05 & edgeR_logFC > 0, "F1-high", fifelse(edgeR_FDR < 0.05 & edgeR_logFC < 0, "Other", "NS"))]
labels <- unique(c(head(volcano[edgeR_logFC > 0][order(edgeR_PValue)]$gene, 10), head(volcano[edgeR_logFC < 0][order(edgeR_PValue)]$gene, 10)))
p2 <- ggplot(volcano, aes(edgeR_logFC, neglog10FDR, colour = class)) +
  geom_point(size = 0.7, alpha = 0.65) +
  geom_text(data = volcano[gene %in% labels], aes(label = gene), size = 2.3, check_overlap = TRUE, vjust = -0.5, show.legend = FALSE) +
  scale_colour_manual(values = c("F1-high" = "#C94755", "Other" = "#3C78A8", "NS" = "#B8B8B8")) +
  labs(x = "edgeR log2 fold change", y = "-log10 FDR", title = "F1-high primary CNV clone vs other primary clones") +
  theme_classic(base_size = 9)
ggsave(file.path(out, "figures", "F1_high_clone_pseudobulk_DE_volcano_v1.png"), p2, width = 6.5, height = 4.8, dpi = 300, bg = "white")
ggsave(file.path(out, "figures", "F1_high_clone_pseudobulk_DE_volcano_v1.pdf"), p2, width = 6.5, height = 4.8)

heat_genes <- unique(c(head(de[edgeR_logFC > 0][order(edgeR_PValue)]$gene, 15), head(de[edgeR_logFC < 0][order(edgeR_PValue)]$gene, 15)))
heat <- as.data.table(deltas[heat_genes, , drop = FALSE], keep.rownames = "gene")
heat <- melt(heat, id.vars = "gene", variable.name = "curated_patient_id", value.name = "logCPM_delta")
heat[, gene := factor(gene, levels = rev(heat_genes))]
p3 <- ggplot(heat, aes(curated_patient_id, gene, fill = logCPM_delta)) +
  geom_tile() + scale_fill_gradient2(low = "#355F9D", mid = "white", high = "#C94755", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "logCPM delta", title = "Patient-paired expression differences") +
  theme_minimal(base_size = 8) + theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, size = 6.5))
ggsave(file.path(out, "figures", "F1_high_clone_top_DE_patient_heatmap_v1.png"), p3, width = 8.2, height = 6.8, dpi = 300, bg = "white")
ggsave(file.path(out, "figures", "F1_high_clone_top_DE_patient_heatmap_v1.pdf"), p3, width = 8.2, height = 6.8)

metadata <- list(
  schema_version = "part1_cnv_clone_cnmf_function_v1",
  completed_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  n_input_cells = length(cells),
  n_input_genes = length(genes),
  n_evaluable_patients = length(patient_levels),
  n_pseudobulk_samples = nrow(sample_meta),
  n_edgeR_tested_genes = nrow(de),
  n_rank7_programs = 7L,
  pathway_sources = c("MSigDB Hallmark v5.2 bundled with VISION", "scMetabolism KEGG_metabolism_nc.gmt", "scMetabolism REACTOME_metabolism.gmt"),
  statistical_unit = "patient-blocked paired pseudobulk",
  interpretation_limits = c("F1 was preselected with TYMP evidence", "CNV clones are expression-inferred and not mutation-validated", "manual figure review is recorded separately"),
  dependency_changes = FALSE
)
jsonlite::write_json(metadata, file.path(out, "run_metadata_v1.json"), auto_unbox = TRUE, pretty = TRUE)
writeLines(capture.output(sessionInfo()), file.path(out, "session_info_v1.txt"))

manifest_files <- setdiff(list.files(out, recursive = TRUE, full.names = TRUE), c(file.path(out, "RUN_COMPLETE"), file.path(out, "file_manifest_v1.tsv")))
manifest_files <- manifest_files[file.info(manifest_files)$isdir == FALSE]
sha256_file <- function(path) strsplit(system2("sha256sum", path, stdout = TRUE), "[[:space:]]+")[[1]][1]
manifest <- data.table(
  file = substring(manifest_files, nchar(out) + 2),
  size_bytes = file.info(manifest_files)$size,
  sha256 = vapply(manifest_files, sha256_file, character(1))
)
fwrite(manifest, file.path(out, "file_manifest_v1.tsv"), sep = "\t")
writeLines(format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"), file.path(out, "RUN_COMPLETE"))
cat("PART1_CNV_CLONE_CNMF_FUNCTION_V1_PASS\n")
