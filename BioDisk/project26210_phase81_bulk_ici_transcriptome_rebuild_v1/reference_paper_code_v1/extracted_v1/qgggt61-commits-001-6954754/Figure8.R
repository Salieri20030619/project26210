suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
infile <- "ELisa.xlsx"                 
outdir <- "mini_heatmap_out_CM769"
dir.create(outdir, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})
outdir <- "mini_heatmap_out_CM_fixed3"
dir.create(outdir, showWarnings = FALSE)
sheets_to_plot <- c("769_Elisa", "786_Elisa")
tile_mm  <- 4
clip_lim <- c(-2, 2)
map_769 <- list(
  "CM-769"         = c("CM-769_1","CM-769_2","CM-769_3"),
  "CM-769+TGFB1"   = c("CM-769+TGFB1_1","CM-769+TGFB1_2","CM-769+TGFB1_3"),
  "CM-769+POSTNab" = c("CM-769+POSTNab_1","CM-769+POSTNab_2","CM-769+POSTNab_3")
)
map_786 <- list(
  "CM-786"         = c("CM-786_1","CM-786_2","CM-786_3"),
  "CM-786+TGFB1"   = c("CM-786+TGFB1_1","CM-786+TGFB1_2","CM-786+TGFB1_3"),
  "CM-786+POSTNab" = c("CM-786+POSTNab_1","CM-786+POSTNab_2","CM-786+POSTNab_3")
)

  
  Heatmap(
    mat_z,
    name = "Z",
    col  = col_fun,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    column_split    = grp,         #
    top_annotation  = top_anno,
    width  = unit(max(1, ncol(mat_z)) * tile_mm, "mm"),
    height = unit(max(1, nrow(mat_z)) * tile_mm, "mm"),
    rect_gp = gpar(col = "white", lwd = 0.8),
    show_row_names     = TRUE,
    row_names_gp       = gpar(fontsize = 9, fontface = "italic"),
    show_column_names  = TRUE,
    column_labels      = col_labels,
    column_names_rot   = 0,
    column_names_gp    = gpar(fontsize = 8),
    border = FALSE,
    column_title    = sheet_name,
    column_title_gp = gpar(fontsize = 11, fontface = "bold")
  )
}

if (length(ht_list) >= 2) {
  ht_combo <- Reduce(`+`, ht_list[1:2])
  pdf(file.path(outdir, "heatmap_combo_first2.pdf"), width = 12, height = 4)
  draw(ht_combo, merge_legend = TRUE)
  dev.off()
  
  png(file.path(outdir, "heatmap_combo_first2.png"), width = 3600, height = 1200, res = 300)
  draw(ht_combo, merge_legend = TRUE)
  dev.off()
}


suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)   # stat_compare_means
  library(ggGantt)  # geom_gantt
})

read_sheet_as_long <- function(file, sheet){
  df <- read_excel(file, sheet = sheet)
  stopifnot(ncol(df) >= 2)
  gene_col <- names(df)[1]
  
  df_long <- df |>
    rename(Gene = !!gene_col) |>
    pivot_longer(-Gene, names_to = "Group_raw", values_to = "Value") |>
    mutate(
      Treatment = str_squish(Treatment),
      Gene      = as.factor(Gene)
    ) |>
    filter(!is.na(Value))
  df_long
}
plot_facets_by_treatment <- function(df_long, sheet_name)
  if (grepl("786", sheet_name, ignore.case = TRUE)){
    ord <- c("CM-786", "CM-786+TGFB1", "CM-786+POSTNab")
  } else if (grepl("769", sheet_name, ignore.case = TRUE)){
    ord <- c("CM-769", "CM-769+TGFB1", "CM-769+POSTNab")
  } else {
    ord <- c("EV","OE","EV+TGFB1","OE+TGFB1ab")
  }
  
  df_long <- df_long |>
    mutate(Treatment = factor(Treatment, levels = ord)) |>
    filter(Treatment %in% ord)
  pal3 <- c(
    "CM-769"           = "#C2A5A1",
    "CM-769+TGFB1"     = "#A8CBB7",
    "CM-769+POSTNab"   = "#A5B8C2",
    "CM-786"           = "#C2A5A1",
    "CM-786+TGFB1"     = "#A8CBB7",
    "CM-786+POSTNab"   = "#A5B8C2",
    "EV"               = "#C2A5A1",
    "OE"               = "#A8CBB7",
    "EV+TGFB1"         = "#A5B8C2",
    "OE+TGFB1ab"       = "#D1C7B7"
  )
  pal_used <- pal3[levels(df_long$Treatment)]
  pal_used[is.na(pal_used)] <- "#bbbbbb"

  w <- max(4, 1.6 * length(unique(df_long$Gene)))
  h <- 8
  
  p <- ggplot(df_long, aes(x = Treatment, y = Value)) +
    geom_point(alpha = 0) +
    geom_gantt(
      aes(fill = Treatment),
      color    = "black",
      stroke   = 0.8,
      width    = 0.5,
      add_line = TRUE,
      position = position_dodge(0.5)
    ) +
    facet_wrap(~ Gene, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = pal_used, drop = FALSE) +
    stat_compare_means(
      method  = "kruskal.test",
      label   = "p.format",
      label.y = y_top
    ) +
    theme_classic(base_size = 14) +
    ylab("Expression") + xlab("") +
    theme(
   
      axis.line.x = element_line(linewidth = 0.1),
      axis.line.y = element_line(linewidth = 0.1),
     
      axis.ticks.x = element_line(linewidth = 0.1),
      axis.ticks.y = element_line(linewidth = 0.1),
      axis.ticks.length = unit(4, "pt"),
  
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      legend.position = "none",
      strip.background = element_blank(),
      strip.text = element_text(face = "bold")
    )
  
  out_pdf <- paste0(sheet_name, ".pdf")
  ggsave(out_pdf, p, width = w, height = h)
  message("Saved: ", normalizePath(out_pdf))
  invisible(p)
}

