#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

input_rds <- Sys.getenv(
  "INPUT_RDS",
  "object/11_malignant_epithelial_cNMF_inferCNV_copyKAT_integrated.rds"
)
out_dir <- Sys.getenv(
  "OUT_DIR",
  "results/scDRS_Factor_TYMP_projection_20260713"
)
score_file <- Sys.getenv("SCDRS_SCORE_FILE", "")
score_name <- Sys.getenv("SCDRS_SCORE_NAME", "scDRS Disease Score")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "rds"), recursive = TRUE, showWarnings = FALSE)

message("[scDRS] Reading malignant epithelial object: ", input_rds)
obj <- readRDS(input_rds)
if (!"Factor_TYMP_score" %in% colnames(obj@meta.data)) {
  stop("Missing Factor_TYMP_score in metadata: ", input_rds)
}
if (!"X_umap" %in% Reductions(obj)) {
  stop("Missing X_umap reduction in object: ", input_rds)
}

get_layer_data <- function(object, assay, layer) {
  out <- tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e) NULL
  )
  if (is.null(out)) {
    out <- tryCatch(GetAssayData(object, assay = assay, slot = layer), error = function(e) NULL)
  }
  out
}

zscore <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
  as.numeric(scale(x))
}

read_external_scdrs <- function(path, cells) {
  dt <- fread(path)
  cell_candidates <- c("cell", "cell_id", "barcode", "Cell", "CELL", "index")
  score_candidates <- c(
    "scDRS_Disease_Score", "scdrs_disease_score", "scDRS_score",
    "scdrs_score", "disease_score", "Disease_Score", "score", "zscore"
  )
  cell_col <- intersect(cell_candidates, names(dt))[1]
  if (is.na(cell_col)) {
    cell_col <- names(dt)[1]
  }
  score_col <- intersect(score_candidates, names(dt))[1]
  if (is.na(score_col)) {
    numeric_cols <- names(dt)[vapply(dt, is.numeric, logical(1))]
    numeric_cols <- setdiff(numeric_cols, cell_col)
    if (!length(numeric_cols)) {
      stop("No numeric scDRS score column found in: ", path)
    }
    score_col <- numeric_cols[1]
  }
  out <- dt[, .(
    cell = as.character(get(cell_col)),
    scDRS_Disease_Score = as.numeric(get(score_col))
  )]
  out <- out[cell %in% cells]
  out
}

make_program_score <- function(object) {
  assay <- if ("RNA" %in% Assays(object)) "RNA" else DefaultAssay(object)
  expr <- get_layer_data(object, assay, "data")
  if (is.null(expr) || nrow(expr) == 0) {
    message("[scDRS] No normalized data layer found; running LogNormalize on assay ", assay)
    DefaultAssay(object) <- assay
    object <- NormalizeData(object, assay = assay, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    expr <- get_layer_data(object, assay, "data")
  }
  genes <- unique(c(
    "VHL", "MET", "EPAS1", "EGLN1", "EGLN2", "EGLN3",
    "FH", "FLCN", "BAP1", "PBRM1", "SETD2", "KDM5C", "KDM6A",
    "MTOR", "TSC1", "TSC2", "PTEN", "SDHA", "SDHB", "SDHC", "SDHD",
    "MITF", "CA9", "NDUFA4L2", "SLC2A1", "VEGFA", "PDK1", "LDHA",
    "ALDOA", "ENO1", "PGK1", "GAPDH", "TPI1", "HK2"
  ))
  genes <- intersect(genes, rownames(expr))
  if (length(genes) < 5) {
    stop("Too few RCC/scDRS-like genes found in expression matrix: ", length(genes))
  }
  raw_score <- Matrix::colMeans(expr[genes, colnames(object), drop = FALSE])
  list(
    object = object,
    score = as.numeric(raw_score),
    genes = genes
  )
}

cells <- colnames(obj)
score_source <- "external_scDRS"
used_genes <- character()

if (nzchar(score_file) && file.exists(score_file)) {
  message("[scDRS] Reading external scDRS score table: ", score_file)
  sc <- read_external_scdrs(score_file, cells)
  meta_score <- data.table(cell = cells)
  meta_score <- merge(meta_score, sc, by = "cell", all.x = TRUE, sort = FALSE)
  if (sum(!is.na(meta_score$scDRS_Disease_Score)) < 100) {
    stop("External scDRS table matched too few cells: ", sum(!is.na(meta_score$scDRS_Disease_Score)))
  }
  score <- meta_score$scDRS_Disease_Score
  score_method <- paste0("External scDRS table: ", basename(score_file))
} else {
  message("[scDRS] No external scDRS score table found; computing provisional RCC GWAS-like program score.")
  program <- make_program_score(obj)
  obj <- program$object
  score <- program$score
  used_genes <- program$genes
  score_source <- "provisional_RCC_GWAS_like_program"
  score_method <- paste0(
    "Provisional RCC GWAS-like disease program score from ",
    length(used_genes),
    " RCC susceptibility / ccRCC malignant-program genes; TYMP excluded."
  )
}

obj$scDRS_Disease_Score <- score
obj$scDRS_Disease_Score_z <- zscore(score)
obj$scDRS_score_source <- score_source

umap <- as.data.table(Embeddings(obj, "X_umap"), keep.rownames = "cell")
setnames(umap, names(umap)[2:3], c("UMAP_1", "UMAP_2"))
meta <- as.data.table(obj@meta.data, keep.rownames = "cell")
plot_dt <- merge(
  umap,
  meta[, .(
    cell,
    Factor_TYMP_score,
    Factor_TYMP_score_z,
    scDRS_Disease_Score,
    scDRS_Disease_Score_z,
    dominant_factor,
    patient_id,
    sample_id,
    copykat_cna_burden,
    highest_malignancy_clone
  )],
  by = "cell",
  all.x = TRUE,
  sort = FALSE
)
plot_dt <- plot_dt[is.finite(scDRS_Disease_Score_z) & is.finite(Factor_TYMP_score_z)]

ct <- cor.test(plot_dt$Factor_TYMP_score_z, plot_dt$scDRS_Disease_Score_z, method = "pearson")
r_lab <- sprintf("Pearson r = %.3f", unname(ct$estimate))
p_lab <- if (ct$p.value < 2.2e-16) "P < 2.2e-16" else paste0("P = ", signif(ct$p.value, 3))
display_score_name <- if (score_source == "external_scDRS") {
  score_name
} else {
  "RCC disease risk score"
}
display_score_subtitle <- if (score_source == "external_scDRS") {
  "Malignant epithelial UMAP"
} else {
  "Provisional score"
}
display_score_name_short <- if (score_source == "external_scDRS") {
  "scDRS\nscore (z)"
} else {
  "RCC risk\nscore (z)"
}

pal_score <- c("#05010A", "#27115F", "#7B2D8E", "#D34D73", "#F6A15A", "#FFF4A3")
theme_projection <- function(base_size = 10) {
  theme_classic(base_size = base_size) +
    theme(
      axis.line = element_line(linewidth = 0.35),
      axis.ticks = element_line(linewidth = 0.3),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey30"),
      legend.title = element_text(face = "bold"),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.25),
      panel.grid.minor = element_blank()
    )
}

score_limits <- quantile(plot_dt$scDRS_Disease_Score_z, c(0.01, 0.99), na.rm = TRUE)
p_umap <- ggplot(plot_dt, aes(UMAP_1, UMAP_2, color = scDRS_Disease_Score_z)) +
  geom_point(size = 0.42, alpha = 0.88, stroke = 0) +
  scale_color_gradientn(
    colors = pal_score,
    limits = score_limits,
    oob = squish,
    name = display_score_name_short
  ) +
  coord_equal() +
  labs(
    title = display_score_name,
    subtitle = display_score_subtitle,
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme_void(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, color = "grey35"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(8, 8, 8, 12),
    plot.background = element_rect(fill = "white", color = NA)
  )

p_scatter <- ggplot(plot_dt, aes(Factor_TYMP_score_z, scDRS_Disease_Score_z)) +
  geom_point(aes(color = copykat_cna_burden), size = 0.55, alpha = 0.45, stroke = 0) +
  geom_smooth(method = "lm", se = TRUE, color = "#111111", fill = "#BDBDBD", linewidth = 0.65) +
  scale_color_gradientn(colors = c("#2C7FB8", "#F7F7F7", "#B2182B"), name = "CNA burden") +
  annotate(
    "label",
    x = quantile(plot_dt$Factor_TYMP_score_z, 0.04, na.rm = TRUE),
    y = quantile(plot_dt$scDRS_Disease_Score_z, 0.96, na.rm = TRUE),
    hjust = 0,
    label = paste(r_lab, p_lab, sep = "\n"),
    linewidth = 0.25,
    fill = "white",
    color = "#111111",
    size = 3.4
  ) +
  labs(
    title = "Factor-TYMP links to disease risk score",
    x = "Factor-TYMP score (z)",
    y = paste0(display_score_name, " (z)")
  ) +
  theme_projection(10) +
  theme(plot.margin = margin(8, 12, 8, 8))

p_combo <- plot_grid(p_umap, p_scatter, nrow = 1, rel_widths = c(1, 1.05))

save_plot <- function(plot, stem, width, height) {
  ggsave(file.path(out_dir, "figures", paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf, bg = "white")
  ggsave(file.path(out_dir, "figures", paste0(stem, ".png")), plot, width = width, height = height, dpi = 450, bg = "white")
}

save_plot(p_umap, "A_scDRS_score_UMAP_projection", 4.8, 4.2)
save_plot(p_scatter, "B_Factor_TYMP_vs_scDRS_score_scatter", 4.8, 4.2)
save_plot(p_combo, "Figure_2_3_scDRS_Factor_TYMP_projection", 10.4, 4.5)

fwrite(plot_dt, file.path(out_dir, "tables", "cell_Factor_TYMP_scDRS_scores.tsv"), sep = "\t")
summary_dt <- data.table(
  n_cells = nrow(plot_dt),
  score_source = score_source,
  score_method = score_method,
  n_score_genes = length(used_genes),
  score_genes = paste(used_genes, collapse = ","),
  pearson_r = unname(ct$estimate),
  pearson_p = ct$p.value,
  factor_tymp_score_mean = mean(plot_dt$Factor_TYMP_score_z, na.rm = TRUE),
  scdrs_score_mean = mean(plot_dt$scDRS_Disease_Score_z, na.rm = TRUE)
)
fwrite(summary_dt, file.path(out_dir, "tables", "scDRS_Factor_TYMP_projection_summary.tsv"), sep = "\t")
saveRDS(obj, file.path(out_dir, "rds", "malignant_epithelial_with_scDRS_projection_scores.rds"))

message("[scDRS] Done.")
message("[scDRS] Figures: ", file.path(out_dir, "figures"))
message("[scDRS] Pearson: ", r_lab, "; ", p_lab)
message("[scDRS] Score source: ", score_source)
