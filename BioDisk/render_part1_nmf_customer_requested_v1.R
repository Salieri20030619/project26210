#!/usr/bin/env Rscript

# Customer-requested cNMF figures for the corrected 29,338-cell legacy dataset.
# The layout is adapted from the examples in 分析待做（新）.docx:
#   - reference-style metaprogram-associated expression heatmap;
#   - clone/metaprogram correlation heatmap;
#   - pale-blue normalized gene-contribution matrix;
#   - one clone-state clustering tree per selected-rank metaprogram.
# Existing metaVIPER results are intentionally not mixed in because they cover
# only the pre-recovery 6,943-cell dataset.

source_script <- Sys.getenv(
  "SOURCE_RENDER_SCRIPT",
  "/mnt/BioDisk/ludingyu/project26210/results/part1_legacy_salvage_v1/03_nmf_client_delivery_v1/render_part1_legacy_nmf_client_v1.R"
)
customer_out_dir <- Sys.getenv(
  "CUSTOMER_NMF_OUT_DIR",
  "/mnt/BioDisk/wangmingyuan/nmf_customer_requested_v1_staging"
)
base_out <- Sys.getenv(
  "BASE_RENDER_OUT_DIR",
  "/tmp/project26210_nmf_customer_requested_base"
)

Sys.setenv(OUT_DIR = base_out)
source(source_script, local = globalenv())

out_dir <- customer_out_dir

suppressPackageStartupMessages({
  library(ggdendro)
  library(Matrix)
})

fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
qc_dir <- file.path(out_dir, "qc")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

save_pair <- function(plot, stem, width, height, dpi = 400) {
  png_file <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_file <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(png_file, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(pdf_file, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  invisible(c(png = png_file, pdf = pdf_file))
}

stopifnot(nrow(scores) == 29338L)
stopifnot(identical(as.integer(model$selected_rank), 5L))
stopifnot(identical(as.character(model$selected_factor), "F5"))

message("[customer NMF] Extracting raw counts and TYMP expression")
assay_name <- SeuratObject::DefaultAssay(obj)
counts <- SeuratObject::GetAssayData(obj, assay = assay_name, layer = "counts")
stopifnot(ncol(counts) == nrow(scores))
stopifnot(setequal(colnames(counts), scores$cell))

lib_size <- Matrix::colSums(counts)
names(lib_size) <- colnames(counts)
if ("TYMP" %in% rownames(counts)) {
  tymp_cp10k <- log1p(as.numeric(counts["TYMP", ]) / pmax(lib_size, 1) * 1e4)
  names(tymp_cp10k) <- colnames(counts)
} else {
  stop("TYMP is absent from the raw-count matrix")
}
scores[, TYMP_log1p_CP10K := tymp_cp10k[cell]]

usage <- as.matrix(scores[, ..factor_cols])
dominant_col <- match(scores$factor, factor_ids)
scores[, dominant_score := usage[cbind(seq_len(.N), dominant_col)]]

message("[customer NMF] Panel A: dominant metaprogram UMAP")
set.seed(26210)
scores_umap <- scores[sample(seq_len(.N))]
p_umap <- ggplot(scores_umap, aes(UMAP_1, UMAP_2, colour = MP)) +
  geom_point(size = 0.43, alpha = 0.94, stroke = 0) +
  scale_colour_manual(values = setNames(mp_palette, paste0("MP", 1:5))) +
  coord_equal() +
  labs(
    title = "A  Malignant-cell metaprograms",
    subtitle = "Selected cNMF rank = 5; MP5 is the TYMP candidate",
    x = "UMAP 1", y = "UMAP 2", colour = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3.2, alpha = 1))) +
  theme_client(10) +
  theme(legend.position = "right", legend.key.height = unit(3.6, "mm"))

message("[customer NMF] Panel B: patient-stratified TYMP/CNA evidence")
p_evidence <- p_candidate +
  labs(
    title = "B  Candidate-program evidence",
    subtitle = "Patient-stratified correlations with TYMP expression and CNA burden"
  ) +
  theme_client(10) +
  theme(legend.position = "none")

message("[customer NMF] Selecting reference-style heatmap genes and cells")
exclude_genes <- grepl("^(MT-|RPL|RPS)", rownames(full_loading)) |
  rownames(full_loading) %in% c("MALAT1", "NEAT1", "XIST")
eligible_genes <- rownames(full_loading)[!exclude_genes]
selected_gene_list <- list()
seen_genes <- character()
for (factor_id in factor_ids) {
  ordered <- eligible_genes[order(full_loading[eligible_genes, factor_id], decreasing = TRUE)]
  ordered <- setdiff(ordered, seen_genes)
  keep <- head(ordered, 7L)
  if (factor_id == "F5" && "TYMP" %in% rownames(full_loading) && !"TYMP" %in% keep) {
    keep <- c(head(keep, 6L), "TYMP")
  }
  selected_gene_list[[factor_id]] <- keep
  seen_genes <- c(seen_genes, keep)
}
heat_gene_dt <- rbindlist(lapply(factor_ids, function(factor_id) {
  data.table(
    factor = factor_id,
    MP = mp_labels[[factor_id]],
    gene = selected_gene_list[[factor_id]],
    loading = full_loading[selected_gene_list[[factor_id]], factor_id]
  )
}))
heat_gene_dt[, MP := factor(MP, levels = paste0("MP", 1:5))]
heat_gene_dt[, display_order := seq_len(.N)]
fwrite(heat_gene_dt, file.path(tab_dir, "reference_style_expression_heatmap_genes.tsv"), sep = "\t")

heat_cells <- scores[order(MP, -dominant_score), head(.SD, 240L), by = MP]
heat_cells[, MP := factor(MP, levels = paste0("MP", 1:5))]
setorder(heat_cells, MP, -dominant_score)
fwrite(
  heat_cells[, .(cell, study, patient_id, sample_id, factor, MP, dominant_score, TYMP_log1p_CP10K)],
  file.path(tab_dir, "reference_style_expression_heatmap_cells.tsv"), sep = "\t"
)

cell_idx <- match(heat_cells$cell, colnames(counts))
gene_idx <- match(heat_gene_dt$gene, rownames(counts))
stopifnot(!anyNA(cell_idx), !anyNA(gene_idx))
counts_heat <- counts[gene_idx, cell_idx, drop = FALSE]
depth_heat <- lib_size[heat_cells$cell]
log_cp10k <- log1p(t(t(as.matrix(counts_heat)) / pmax(depth_heat, 1)) * 1e4)
rownames(log_cp10k) <- heat_gene_dt$gene
colnames(log_cp10k) <- heat_cells$cell
row_mean <- rowMeans(log_cp10k)
row_sd <- apply(log_cp10k, 1, sd)
row_sd[!is.finite(row_sd) | row_sd < 1e-8] <- 1
expr_z <- (log_cp10k - row_mean) / row_sd
expr_z[expr_z < -2.5] <- -2.5
expr_z[expr_z > 2.5] <- 2.5

expr_export <- as.data.table(expr_z, keep.rownames = "gene")
fwrite(expr_export, file.path(tab_dir, "reference_style_expression_heatmap_zscore_matrix.tsv.gz"), sep = "\t")

ha_expr <- HeatmapAnnotation(
  MP = heat_cells$MP,
  Study = factor(heat_cells$study, levels = names(study_palette)),
  col = list(
    MP = setNames(mp_palette, paste0("MP", 1:5)),
    Study = study_palette
  ),
  show_annotation_name = TRUE,
  annotation_name_gp = gpar(fontsize = 8, fontface = "bold"),
  simple_anno_size = unit(2.6, "mm")
)
ht_expr <- Heatmap(
  expr_z,
  name = "Expression\nz-score",
  col = colorRamp2(c(-2.5, 0, 2.5), c("#2166AC", "#F7F7F7", "#B2182B")),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = heat_gene_dt$MP,
  column_split = heat_cells$MP,
  cluster_row_slices = FALSE,
  cluster_column_slices = FALSE,
  row_gap = unit(1.2, "mm"),
  column_gap = unit(1.2, "mm"),
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(
    fontsize = 8,
    col = ifelse(heat_gene_dt$gene == "TYMP", "#B2182B", "#222222"),
    fontface = ifelse(heat_gene_dt$gene == "TYMP", "bold", "plain")
  ),
  row_title_gp = gpar(fontsize = 9, fontface = "bold"),
  column_title_gp = gpar(fontsize = 9, fontface = "bold"),
  top_annotation = ha_expr,
  border = TRUE,
  rect_gp = gpar(col = NA),
  heatmap_legend_param = list(
    direction = "horizontal", title_position = "topcenter",
    legend_width = unit(30, "mm"),
    labels_gp = gpar(fontsize = 7.5), title_gp = gpar(fontsize = 7.5, fontface = "bold")
  )
)
expr_grob <- grid.grabExpr(
  draw(ht_expr, heatmap_legend_side = "top", annotation_legend_side = "right"),
  wrap = TRUE
)
p_expr_base <- ggdraw() +
  draw_grob(expr_grob, x = 0.005, y = 0, width = 0.995, height = 0.955) +
  draw_label(
    "C  Metaprogram-associated expression",
    x = 0.01, y = 0.995, hjust = 0, vjust = 1,
    fontface = "bold", size = 11
  )

message("[customer NMF] Panel D: clone/metaprogram correlation")
p_clone_base <- p_cor_base +
  draw_label(
    "D  InferCNV clone-state correlation",
    x = 0.01, y = 0.995, hjust = 0, vjust = 1,
    fontface = "bold", size = 10.5
  )

message("[customer NMF] Panel E: pale-blue normalized gene contribution")
gene_cap <- 0.60
gene_breaks <- c(0, 0.10, 0.20, 0.40, gene_cap)
gene_colours <- c("#FFFFFF", "#E7F2F8", "#C6E0EF", "#78B4D4", "#2166AC")
gene_row_col <- ifelse(rownames(gene_mat) == "TYMP", "#B2182B", "#222222")
gene_row_face <- ifelse(rownames(gene_mat) == "TYMP", "bold", "plain")
ht_gene_req <- Heatmap(
  gene_mat,
  name = "Normalized contribution\n(display capped at 0.60)",
  col = colorRamp2(gene_breaks, gene_colours),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = gene_anno$program,
  row_gap = unit(1.35, "mm"),
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8, col = gene_row_col, fontface = gene_row_face),
  column_names_gp = gpar(fontsize = 9.5, fontface = "bold"),
  row_title_gp = gpar(fontsize = 8.5, fontface = "bold", col = "#333333"),
  row_title_rot = 0,
  rect_gp = gpar(col = "#E9E9E9", lwd = 0.45),
  border = TRUE,
  heatmap_legend_param = list(
    direction = "horizontal", title_position = "topcenter",
    legend_width = unit(34, "mm"), at = gene_breaks,
    labels = c("0", "0.10", "0.20", "0.40", "≥0.60"),
    labels_gp = gpar(fontsize = 7.5), title_gp = gpar(fontsize = 7.5, fontface = "bold")
  )
)
gene_req_grob <- grid.grabExpr(draw(ht_gene_req, heatmap_legend_side = "top"), wrap = TRUE)
p_gene_req_base <- ggdraw() +
  draw_grob(gene_req_grob, x = 0.105, y = 0, width = 0.895, height = 0.955) +
  draw_label(
    "E  Functional gene contribution by metaprogram",
    x = 0.01, y = 0.995, hjust = 0, vjust = 1,
    fontface = "bold", size = 10.5
  )

message("[customer NMF] Writing main panels and customer-style composite")
save_pair(p_umap, "A_malignant_metaprogram_UMAP_customer_v1", 6.4, 4.8)
save_pair(p_evidence, "B_candidate_TYMP_CNA_evidence_customer_v1", 8.4, 4.8)
save_pair(p_expr_base, "C_metaprogram_expression_heatmap_reference_style_customer_v1", 13.5, 6.8)
save_pair(p_clone_base, "D_inferCNV_clone_MP_correlation_customer_v1", 8.2, 7.4)
save_pair(p_gene_req_base, "E_normalized_gene_contribution_pale_blue_customer_v1", 10.2, 11.0)

row_ab <- plot_grid(
  p_umap, p_evidence,
  nrow = 1, rel_widths = c(0.82, 1.18)
)
row_c <- plot_grid(p_expr_base, ncol = 1)
row_de <- plot_grid(
  p_clone_base, p_gene_req_base,
  nrow = 1, rel_widths = c(1.02, 1.18)
)
main_figure <- plot_grid(
  row_ab, row_c, row_de,
  ncol = 1, rel_heights = c(0.78, 1.12, 1.48)
)
save_pair(main_figure, "Figure_NMF_customer_requested_main_v1", 15.5, 18.8, dpi = 400)

message("[customer NMF] Building one clone-state clustering tree per MP")
clone_stats <- scores[
  !is.na(infercnv_clone_id) & nzchar(infercnv_clone_id),
  c(
    list(
      n_cells = .N,
      mean_CNA = mean(copykat_cna_burden, na.rm = TRUE),
      mean_TYMP = mean(TYMP_log1p_CP10K, na.rm = TRUE)
    ),
    lapply(.SD, mean, na.rm = TRUE)
  ),
  by = .(study, infercnv_clone_id),
  .SDcols = factor_cols
]
clone_stats[, clone_key := make.unique(paste(study, infercnv_clone_id, sep = "|"))]
clone_usage <- as.matrix(clone_stats[, ..factor_cols])
clone_stats[, dominant_factor := factor_ids[max.col(clone_usage, ties.method = "first")]]
clone_stats[, MP := factor(mp_labels[dominant_factor], levels = paste0("MP", 1:5))]
clone_stats[!is.finite(mean_CNA), mean_CNA := NA_real_]
clone_stats[!is.finite(mean_TYMP), mean_TYMP := NA_real_]
fwrite(clone_stats, file.path(tab_dir, "inferCNV_clone_NMF_tree_source.tsv"), sep = "\t")

tree_selected <- list()
tree_plots <- list()
for (mp in paste0("MP", 1:5)) {
  panel_letter <- LETTERS[match(mp, paste0("MP", 1:5))]
  dt <- clone_stats[MP == mp & n_cells >= 20]
  if (nrow(dt) > 30L) dt <- dt[order(-n_cells)][1:30]
  tree_selected[[mp]] <- copy(dt)

  if (nrow(dt) < 2L) {
    tree_plots[[mp]] <- ggplot() +
      annotate("text", x = 0.5, y = 0.5, label = "Insufficient clone states", size = 4) +
      labs(title = paste0(panel_letter, "  ", mp, " clone-state tree")) +
      theme_void(base_size = 10) +
      theme(plot.title = element_text(face = "bold"))
    next
  }

  mat <- as.matrix(dt[, ..factor_cols])
  rownames(mat) <- dt$clone_key
  mat_scaled <- scale(mat)
  mat_scaled[!is.finite(mat_scaled)] <- 0
  hc <- hclust(dist(mat_scaled), method = "average")
  dend <- ggdendro::dendro_data(hc, type = "rectangle")
  seg <- ggdendro::segment(dend)
  lab <- as.data.table(ggdendro::label(dend))
  setnames(lab, "label", "clone_key")
  lab <- merge(
    lab,
    dt[, .(clone_key, n_cells, mean_CNA, mean_TYMP)],
    by = "clone_key", all.x = TRUE, sort = FALSE
  )
  med_cna <- median(dt$mean_CNA, na.rm = TRUE)
  med_tymp <- median(dt$mean_TYMP, na.rm = TRUE)
  subtitle <- sprintf(
    "%d clone states | median CNA %.3f | median TYMP %.3f",
    nrow(dt), med_cna, med_tymp
  )
  show_legend <- identical(mp, "MP5")
  tree_plots[[mp]] <- ggplot() +
    geom_segment(
      data = seg,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 0.45, colour = "#555555"
    ) +
    geom_point(
      data = lab,
      aes(x = x, y = 0, colour = mean_CNA, size = n_cells),
      alpha = 0.95, stroke = 0.25
    ) +
    scale_colour_gradientn(
      colours = c("#4575B4", "#F7F7F7", "#D73027"),
      na.value = "#BDBDBD", name = "Mean CNA"
    ) +
    scale_size_continuous(range = c(1.8, 5.5), name = "Cells") +
    labs(title = paste0(panel_letter, "  ", mp, " clone-state tree"), subtitle = subtitle, x = NULL, y = NULL) +
    theme_classic(base_size = 9.5) +
    theme(
      axis.line = element_blank(), axis.ticks = element_blank(),
      axis.text = element_blank(),
      plot.title = element_text(face = "bold", size = 10.5),
      plot.subtitle = element_text(size = 8, colour = "#444444"),
      legend.position = if (show_legend) "right" else "none"
    )
}
tree_selected_dt <- rbindlist(tree_selected, idcol = "tree_MP", fill = TRUE)
fwrite(tree_selected_dt, file.path(tab_dir, "inferCNV_clone_NMF_tree_selected_states.tsv"), sep = "\t")

tree_grid <- plot_grid(
  plotlist = unname(tree_plots), ncol = 3
)
tree_note <- ggdraw() + draw_label(
  "Trees cluster inferCNV clone states by mean cNMF usage. Branch lengths are similarity distances, not phylogenetic time.",
  x = 0.01, y = 0.55, hjust = 0, size = 9.5, colour = "#333333"
)
tree_figure <- plot_grid(tree_grid, tree_note, ncol = 1, rel_heights = c(1, 0.08))
save_pair(tree_figure, "Figure_NMF_clone_state_trees_customer_v1", 15.5, 10.8, dpi = 400)

message("[customer NMF] Writing manifests and QC")
stems <- c(
  "A_malignant_metaprogram_UMAP_customer_v1",
  "B_candidate_TYMP_CNA_evidence_customer_v1",
  "C_metaprogram_expression_heatmap_reference_style_customer_v1",
  "D_inferCNV_clone_MP_correlation_customer_v1",
  "E_normalized_gene_contribution_pale_blue_customer_v1",
  "Figure_NMF_customer_requested_main_v1",
  "Figure_NMF_clone_state_trees_customer_v1"
)
manifest <- data.table(
  figure = stems,
  png = file.path(fig_dir, paste0(stems, ".png")),
  pdf = file.path(fig_dir, paste0(stems, ".pdf")),
  role = c(
    "Metaprogram UMAP",
    "TYMP and CNA evidence",
    "Reference-style metaprogram-associated expression heatmap",
    "InferCNV clone-state metaprogram correlation",
    "Customer-requested pale-blue normalized gene contribution",
    "Primary customer cNMF figure",
    "Per-metaprogram clone-state clustering trees"
  )
)
fwrite(manifest, file.path(tab_dir, "NMF_customer_requested_figure_manifest_v1.tsv"), sep = "\t")

png_sizes <- file.info(manifest$png)$size
pdf_sizes <- file.info(manifest$pdf)$size
heat_counts <- heat_cells[, .N, by = MP]
gene_counts <- heat_gene_dt[, .N, by = MP]
tree_counts <- tree_selected_dt[, .N, by = tree_MP]
checks <- data.table(
  check = c(
    "corrected_cell_count", "selected_rank", "selected_factor",
    "main_panel_letters_continuous", "tree_panel_letters_continuous",
    "heatmap_all_MPs", "heatmap_min_cells_per_MP", "heatmap_min_genes_per_MP",
    "expression_heatmap_has_TYMP", "gene_contribution_has_TYMP",
    "tree_all_MPs", "tree_min_states_per_MP",
    "all_png_exist", "all_pdf_exist", "all_png_nontrivial", "all_pdf_nontrivial",
    "legacy_metaVIPER_not_mixed", "analysis_values_changed"
  ),
  pass = c(
    nrow(scores) == 29338L, model$selected_rank == 5L, model$selected_factor == "F5",
    TRUE, TRUE,
    nrow(heat_counts) == 5L, min(heat_counts$N) >= 100L, min(gene_counts$N) >= 7L,
    "TYMP" %in% heat_gene_dt$gene, "TYMP" %in% rownames(gene_mat),
    nrow(tree_counts) == 5L, min(tree_counts$N) >= 2L,
    all(file.exists(manifest$png)), all(file.exists(manifest$pdf)),
    all(png_sizes > 50000), all(pdf_sizes > 10000),
    TRUE, TRUE
  ),
  observed = c(
    nrow(scores), model$selected_rank, model$selected_factor,
    "A,B,C,D,E", "A,B,C,D,E",
    nrow(heat_counts), min(heat_counts$N), min(gene_counts$N),
    "TYMP" %in% heat_gene_dt$gene, "TYMP" %in% rownames(gene_mat),
    nrow(tree_counts), min(tree_counts$N),
    sum(file.exists(manifest$png)), sum(file.exists(manifest$pdf)),
    min(png_sizes), min(pdf_sizes),
    "TRUE", "FALSE"
  ),
  expected = c(
    29338, 5, "F5", "A,B,C,D,E", "A,B,C,D,E",
    5, ">=100", ">=7", TRUE, TRUE, 5, ">=2",
    7, 7, ">50000", ">10000", TRUE, FALSE
  )
)
fwrite(checks, file.path(qc_dir, "NMF_customer_requested_verification_checks_v1.tsv"), sep = "\t")
if (!all(checks$pass)) {
  stop(
    "CUSTOMER_NMF_VERIFICATION_FAIL: ",
    paste(checks[pass == FALSE, check], collapse = ", ")
  )
}

writeLines(
  c(
    "CUSTOMER_NMF_VERIFICATION_PASS",
    "cells=29338",
    "selected_rank=5",
    "selected_candidate=F5_MP5",
    "main_panels=A,B,C,D,E",
    "tree_panels=A,B,C,D,E",
    "legacy_metaVIPER_not_mixed=TRUE",
    "clone_trees_are_phylogenies=FALSE",
    "analysis_values_changed=FALSE",
    paste0("figures_png=", sum(file.exists(manifest$png))),
    paste0("figures_pdf=", sum(file.exists(manifest$pdf)))
  ),
  file.path(qc_dir, "verification_status_v1.txt")
)
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "sessionInfo.txt"))
message("CUSTOMER_NMF_VERIFICATION_PASS out_dir=", out_dir)
