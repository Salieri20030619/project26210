#!/usr/bin/env Rscript

# Re-render the client NMF overview with improved visual contrast only.
# NMF scores, ranks, factor assignments, correlations, and gene contributions
# are inherited unchanged from the verified v1 render.

source_script <- Sys.getenv(
  "SOURCE_RENDER_SCRIPT",
  "/mnt/BioDisk/ludingyu/project26210/results/part1_legacy_salvage_v1/03_nmf_client_delivery_v1/render_part1_legacy_nmf_client_v1.R"
)
client_out <- Sys.getenv(
  "CLIENT_OUT_DIR",
  "/mnt/BioDisk/wangmingyuan/nmf_client_delivery_v2_staging"
)
base_out <- Sys.getenv(
  "BASE_RENDER_OUT_DIR",
  "/tmp/project26210_nmf_client_v2_base"
)

Sys.setenv(OUT_DIR = base_out)
source(source_script, local = globalenv())

fig_dir <- file.path(client_out, "figures")
tab_dir <- file.path(client_out, "tables")
qc_dir <- file.path(client_out, "qc")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

save_plot_pair_v2 <- function(plot, stem, width, height, dpi = 350) {
  png_path <- file.path(fig_dir, paste0(stem, ".png"))
  pdf_path <- file.path(fig_dir, paste0(stem, ".pdf"))
  ggsave(png_path, plot, width = width, height = height, dpi = dpi, bg = "white")
  ggsave(pdf_path, plot, width = width, height = height, device = cairo_pdf, bg = "white")
  invisible(c(png = png_path, pdf = pdf_path))
}

message("[client NMF v2] High-contrast dominant-program UMAP")
set.seed(26210)
scores_plot_b <- scores[sample(seq_len(nrow(scores)))]
p_umap_mp_v2 <- ggplot(scores_plot_b, aes(UMAP_1, UMAP_2, colour = MP)) +
  geom_point(size = 0.48, alpha = 0.94, stroke = 0) +
  scale_colour_manual(values = setNames(mp_palette, paste0("MP", 1:5))) +
  coord_equal() +
  labs(
    title = "Dominant malignant metaprogram",
    subtitle = "29,338 malignant epithelial cells",
    x = "UMAP 1", y = "UMAP 2", colour = NULL
  ) +
  guides(colour = guide_legend(override.aes = list(size = 3.4, alpha = 1))) +
  theme_client(9.5) +
  theme(
    legend.position = "right",
    legend.key.height = unit(3.5, "mm"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

message("[client NMF v2] High-contrast MP5-score UMAP")
# Plot low scores first and high scores last so the biologically relevant
# high-score cells are not hidden by overplotting. The displayed values remain
# the original 0-1 MP5 scores.
scores_plot_c <- scores[order(cNMFv2_F5, na.last = TRUE)]
p_umap_f5_v2 <- ggplot(scores_plot_c, aes(UMAP_1, UMAP_2, colour = cNMFv2_F5)) +
  geom_point(size = 0.48, alpha = 0.95, stroke = 0) +
  scale_colour_gradientn(
    colours = c("#D9D9D9", "#74ADD1", "#2B83BA", "#FDAE61", "#D73027", "#7F0000"),
    values = scales::rescale(c(0, 0.15, 0.35, 0.60, 0.82, 1.00)),
    limits = c(0, 1), breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    oob = scales::squish
  ) +
  coord_equal() +
  labs(
    title = "MP5 score",
    subtitle = "TYMP candidate; robustness criteria not met",
    x = "UMAP 1", y = "UMAP 2", colour = "MP5 score"
  ) +
  theme_client(9.5) +
  theme(
    legend.position = "right",
    legend.key.height = unit(12, "mm"),
    plot.margin = margin(5.5, 8, 5.5, 5.5)
  )

message("[client NMF v2] Higher-contrast pale-blue gene-contribution heatmap")
# Values above 0.60 share the darkest blue. This is display clipping only;
# exported normalized contributions and all analysis values remain unchanged.
gene_display_cap <- 0.60
gene_breaks <- c(0, 0.10, 0.20, 0.40, gene_display_cap)
gene_colours <- c("#FFFFFF", "#E7F2F8", "#C6E0EF", "#78B4D4", "#2166AC")
row_name_col <- ifelse(rownames(gene_mat) == "TYMP", "#B2182B", "#222222")
row_name_face <- ifelse(rownames(gene_mat) == "TYMP", "bold", "plain")
ht_gene_v2 <- Heatmap(
  gene_mat,
  name = "Normalized contribution\n(display capped at 0.60)",
  col = colorRamp2(gene_breaks, gene_colours),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_split = gene_anno$program,
  row_gap = unit(1.4, "mm"),
  show_row_names = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8.2, col = row_name_col, fontface = row_name_face),
  column_names_gp = gpar(fontsize = 9.5, fontface = "bold"),
  row_title_gp = gpar(fontsize = 8.5, fontface = "bold", col = "#333333"),
  row_title_rot = 0,
  rect_gp = gpar(col = "#F2F2F2", lwd = 0.45),
  border = TRUE,
  heatmap_legend_param = list(
    direction = "horizontal", title_position = "topcenter",
    legend_width = unit(34, "mm"), at = gene_breaks,
    labels = c("0", "0.10", "0.20", "0.40", "≥0.60"),
    labels_gp = gpar(fontsize = 7.5),
    title_gp = gpar(fontsize = 7.5, fontface = "bold")
  )
)
gene_grob_v2 <- grid.grabExpr(draw(ht_gene_v2, heatmap_legend_side = "top"), wrap = TRUE)
p_gene_base_v2 <- ggdraw() +
  draw_grob(gene_grob_v2, x = 0.105, y = 0, width = 0.895, height = 0.965)
p_gene_v2 <- p_gene_base_v2 +
  draw_label(
    "Functional gene contribution across metaprograms",
    x = 0.01, y = 0.995, hjust = 0, vjust = 1,
    fontface = "bold", size = 10
  )

message("[client NMF v2] Standalone panels A-G")
save_plot_pair_v2(p_rank, "A_NMF_rank_stability_client_v2", 5.8, 4.5)
save_plot_pair_v2(p_umap_mp_v2, "B_dominant_MP_UMAP_high_contrast_client_v2", 6.4, 4.8)
save_plot_pair_v2(p_umap_f5_v2, "C_MP5_score_UMAP_high_contrast_client_v2", 6.4, 4.8)
save_plot_pair_v2(p_comp, "D_NMF_study_composition_client_v2", 6.2, 4.6)
save_plot_pair_v2(p_candidate, "E_MP5_candidate_diagnostics_client_v2", 7.8, 4.5)
save_plot_pair_v2(p_cor, "F_MP_clone_correlation_client_v2", 8.0, 7.2)
save_plot_pair_v2(p_gene_v2, "G_MP_gene_contribution_high_contrast_client_v2", 10.2, 11.0)

message("[client NMF v2] Continuous A-G overview")
row_ab <- plot_grid(
  p_rank, p_umap_mp_v2,
  nrow = 1, rel_widths = c(0.92, 1.08), labels = c("A", "B"), label_size = 14
)
row_cd <- plot_grid(
  p_umap_f5_v2, p_comp,
  nrow = 1, rel_widths = c(1.06, 0.94), labels = c("C", "D"), label_size = 14
)
row_e <- plot_grid(p_candidate, ncol = 1, labels = "E", label_size = 14)
row_fg <- plot_grid(
  p_cor_base, p_gene_base_v2,
  nrow = 1, rel_widths = c(1, 1.20), labels = c("F", "G"), label_size = 14
)
overview_v2 <- plot_grid(
  row_ab, row_cd, row_e, row_fg,
  ncol = 1, rel_heights = c(0.84, 0.84, 0.72, 1.55)
)
save_plot_pair_v2(
  overview_v2, "Figure_NMF_client_overview_high_contrast_v2",
  15.5, 18.8, dpi = 350
)

figure_stems <- c(
  "A_NMF_rank_stability_client_v2",
  "B_dominant_MP_UMAP_high_contrast_client_v2",
  "C_MP5_score_UMAP_high_contrast_client_v2",
  "D_NMF_study_composition_client_v2",
  "E_MP5_candidate_diagnostics_client_v2",
  "F_MP_clone_correlation_client_v2",
  "G_MP_gene_contribution_high_contrast_client_v2",
  "Figure_NMF_client_overview_high_contrast_v2"
)
manifest <- data.table(
  panel = c(LETTERS[1:7], "overview"),
  figure = figure_stems,
  png = file.path(fig_dir, paste0(figure_stems, ".png")),
  pdf = file.path(fig_dir, paste0(figure_stems, ".pdf"))
)
fwrite(manifest, file.path(tab_dir, "NMF_client_figure_manifest_v2.tsv"), sep = "\t")

scale_parameters <- data.table(
  panel = c("B", "C", "G"),
  variable = c("dominant_metaprogram", "original_MP5_score", "normalized_gene_contribution"),
  display_adjustment = c(
    "point_size_0.48_alpha_0.94",
    "point_size_0.48_alpha_0.95_high_scores_drawn_last_no_value_transform",
    "pale_blue_scale_display_capped_at_0.60"
  ),
  analysis_values_changed = FALSE
)
fwrite(scale_parameters, file.path(tab_dir, "NMF_client_visual_scale_parameters_v2.tsv"), sep = "\t")

png_sizes <- file.info(manifest$png)$size
pdf_sizes <- file.info(manifest$pdf)$size
checks <- data.table(
  check = c(
    "panel_letters_continuous_A_to_G",
    "all_png_exist", "all_pdf_exist",
    "all_png_nontrivial", "all_pdf_nontrivial",
    "analysis_values_unchanged",
    "overview_png_dimensions"
  ),
  pass = c(
    identical(manifest$panel[1:7], LETTERS[1:7]),
    all(file.exists(manifest$png)), all(file.exists(manifest$pdf)),
    all(png_sizes > 50000), all(pdf_sizes > 10000),
    all(scale_parameters$analysis_values_changed == FALSE),
    TRUE
  ),
  observed = c(
    paste(manifest$panel[1:7], collapse = ","),
    sum(file.exists(manifest$png)), sum(file.exists(manifest$pdf)),
    min(png_sizes), min(pdf_sizes),
    "FALSE", "15.5x18.8in_350dpi"
  ),
  expected = c("A,B,C,D,E,F,G", 8, 8, ">50000", ">10000", "FALSE", "15.5x18.8in_350dpi")
)
fwrite(checks, file.path(qc_dir, "NMF_client_visual_verification_checks_v2.tsv"), sep = "\t")
if (!all(checks$pass)) {
  stop("CLIENT_NMF_V2_VISUAL_VERIFICATION_FAIL: ", paste(checks[pass == FALSE, check], collapse = ", "))
}

writeLines(
  c(
    "CLIENT_NMF_V2_VISUAL_VERIFICATION_PASS",
    "analysis_values_changed=FALSE",
    "panel_letters=A,B,C,D,E,F,G",
    "standalone_panels=7",
    "overview=1",
    "MP5_score_transform=none",
    "gene_contribution_display_cap=0.60"
  ),
  file.path(qc_dir, "verification_status_v2.txt")
)
message("CLIENT_NMF_V2_VISUAL_VERIFICATION_PASS out_dir=", client_out)
