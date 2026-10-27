setwd('H:/2025年文章/stMTMC.RCC/Figures/related to figure7')
## ===== 依赖 =====
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

## ===== 路径与输出 =====
infile <- "ELisa.xlsx"                 # ← 换成你的文件名
outdir <- "mini_heatmap_out_CM769"
dir.create(outdir, showWarnings = FALSE)
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
})

## ======= 基本参数 =======
infile <- "ELisa.xlsx"  # ← 改成你的文件名
outdir <- "mini_heatmap_out_CM_fixed3"
dir.create(outdir, showWarnings = FALSE)

# 你要画的 sheet 名（按你的文件来改）
sheets_to_plot <- c("769_Elisa", "786_Elisa")

# 单元格大小 & Z 分数截断
tile_mm  <- 4
clip_lim <- c(-2, 2)

# 热图连续色（莫兰迪：冷灰蓝 → 白 → 暖粉）
col_fun <- circlize::colorRamp2(
  c(-2, 0, 2),
  c("#B7C4C8", "#F3F1ED", "#CDB4BD")
)

# 三组分片配色（莫兰迪）
grp_colors <- c(
  "BASE"    = "#B7C4C8",  # CM-xxx
  "TGFB1"   = "#C9BEDA",  # CM-xxx+TGFB1
  "POSTNab" = "#C8D3C5"   # CM-xxx+POSTNab
)

# 轻度清理列名里的奇怪空格
sanitize_names <- function(x){
  x <- gsub("[\u00A0\u2007\u202F]", " ", x); x <- gsub("\\s+", " ", x); trimws(x)
}

## ======= 明确映射：三组固定到具体列名（不用正则） =======
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

get_map_for_sheet <- function(sn){
  if (grepl("769", sn, ignore.case = TRUE)) return(map_769)
  if (grepl("786", sn, ignore.case = TRUE)) return(map_786)
  stop(sprintf("未识别的 sheet：%s，请在 get_map_for_sheet() 中补充映射。", sn))
}

## ======= 画单个 sheet 的函数（不聚类，固定三组分片） =======
make_ht_for_sheet <- function(sheet_name) {
  # 读表
  df <- read_excel(infile, sheet = sheet_name) |> as.data.frame(check.names = FALSE)
  if (ncol(df) < 2) stop(sprintf("【%s】列数 < 2。", sheet_name))
  
  # 行名 + 数值矩阵
  rn <- as.character(df[[1]])
  matdf <- df[, -1, drop = FALSE]
  colnames(matdf) <- sanitize_names(colnames(matdf))
  
  # 仅保留数值列
  num_cols <- vapply(matdf, is.numeric, logical(1))
  if (!any(num_cols)) stop(sprintf("【%s】没有可用的数值列。", sheet_name))
  mat <- as.matrix(matdf[, num_cols, drop = FALSE])
  rownames(mat) <- rn
  
  # 取该 sheet 的三组映射
  grp_map <- get_map_for_sheet(sheet_name)
  
  # 找实际存在的列（按映射顺序拼接）
  col_order <- unlist(grp_map, use.names = FALSE)
  col_order <- intersect(col_order, colnames(mat))
  if (length(col_order) == 0) stop(sprintf("【%s】三组列名都没在数据里找到。", sheet_name))
  
  # 子集 & 重排
  mat <- mat[, col_order, drop = FALSE]
  
  # 与列对齐的组标签
  grp_vec <- rep(names(grp_map), times = vapply(grp_map, function(v) sum(v %in% col_order), integer(1)))
  grp_vec <- grp_vec[match(col_order, unlist(grp_map, use.names = FALSE))]
  
  # 行清洗
  keep  <- rowSums(is.na(mat)) < ncol(mat)
  keep2 <- rowSums(abs(mat), na.rm = TRUE) > 0
  mat   <- mat[keep & keep2, , drop = FALSE]
  if (nrow(mat) == 0) stop(sprintf("【%s】清洗后无可画行。", sheet_name))
  
  # 行向 Z-score + 截断
  mat_z <- t(scale(t(mat))); mat_z[is.na(mat_z)] <- 0
  mat_z <- pmax(pmin(mat_z, clip_lim[2]), clip_lim[1])
  
  # 列标签（S1/S2/S3）
  repn <- suppressWarnings(as.numeric(sub(".*_(\\d+)$", "\\1", colnames(mat_z))))
  col_labels <- ifelse(is.na(repn), "", paste0("S", repn))
  
  # 三组固定顺序
  grp_levels <- names(grp_map)
  grp <- factor(grp_vec, levels = grp_levels)
  
  # 顶部分组色条
  top_anno <- HeatmapAnnotation(
    Group = grp,
    col = list(Group = setNames(
      c(grp_colors["BASE"], grp_colors["TGFB1"], grp_colors["POSTNab"]),
      grp_levels
    )),
    which = "column",
    annotation_legend_param = list(title = "Group")
  )
  
  # 画图（不聚类 + 分片）
  Heatmap(
    mat_z,
    name = "Z",
    col  = col_fun,
    cluster_rows    = FALSE,
    cluster_columns = FALSE,
    column_split    = grp,         # 三组分片；默认会在每片上方标注组名
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

## ======= 主流程：逐个 sheet 生成并导出 =======
ht_list <- list()
for (sn in sheets_to_plot) {
  message(">>> 绘图：", sn)
  ht <- tryCatch(make_ht_for_sheet(sn),
                 error = function(e){ message("  跳过【", sn, "】：", e$message); NULL })
  if (is.null(ht)) next
  ht_list[[sn]] <- ht
  
  pdf(file.path(outdir, paste0("heatmap_", sn, ".pdf")), width = 6, height = 4)
  draw(ht, merge_legend = TRUE)
  dev.off()
  
  png(file.path(outdir, paste0("heatmap_", sn, ".png")), width = 1800, height = 1200, res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
}

# 可选：把前两个并排导出
if (length(ht_list) >= 2) {
  ht_combo <- Reduce(`+`, ht_list[1:2])
  pdf(file.path(outdir, "heatmap_combo_first2.pdf"), width = 12, height = 4)
  draw(ht_combo, merge_legend = TRUE)
  dev.off()
  
  png(file.path(outdir, "heatmap_combo_first2.png"), width = 3600, height = 1200, res = 300)
  draw(ht_combo, merge_legend = TRUE)
  dev.off()
}

cat("完成。输出目录：", normalizePath(outdir), "\n")




## ========= 依赖 =========
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggpubr)   # stat_compare_means
  library(ggGantt)  # geom_gantt
})

## ========= 参数 =========
infile <- "PCR.xlsx"  # ← 换成你的文件名/路径
## ===== 依赖 =====
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(ggpubr)   # stat_compare_means
  library(ggGantt)  # geom_gantt
})



## ===== 读并“分组”清洗：把每组的 _1/_2/_3 归并为一个 Treatment =====
read_sheet_as_long <- function(file, sheet){
  df <- read_excel(file, sheet = sheet)
  stopifnot(ncol(df) >= 2)
  gene_col <- names(df)[1]
  
  df_long <- df |>
    rename(Gene = !!gene_col) |>
    pivot_longer(-Gene, names_to = "Group_raw", values_to = "Value") |>
    mutate(
      # 去掉末尾的 _1/_2/_3，把重复归入一个“Treatment”组
      Treatment = str_replace(Group_raw, "_\\d+$", ""),
      # 去掉奇怪空格
      Treatment = str_squish(Treatment),
      Gene      = as.factor(Gene)
    ) |>
    filter(!is.na(Value))
  df_long
}

## ===== 绘图函数（按“Treatment”分组；每个 Gene 一个面板）=====
plot_facets_by_treatment <- function(df_long, sheet_name){
  
  # 根据 sheet_name 自动设定 3 组顺序（你也可以手动写死）
  if (grepl("786", sheet_name, ignore.case = TRUE)){
    ord <- c("CM-786", "CM-786+TGFB1", "CM-786+POSTNab")
  } else if (grepl("769", sheet_name, ignore.case = TRUE)){
    ord <- c("CM-769", "CM-769+TGFB1", "CM-769+POSTNab")
  } else {
    # 如果是 EV/OE 体系，就用这套（需要时可切换）
    ord <- c("EV","OE","EV+TGFB1","OE+TGFB1ab")
  }
  
  df_long <- df_long |>
    mutate(Treatment = factor(Treatment, levels = ord)) |>
    filter(Treatment %in% ord)
  
  # 莫兰迪配色（3 组）
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
  
  # 顶部 p 值位置（用全局最大值的 1.15 倍）
  y_top <- 1.15 * max(df_long$Value, na.rm = TRUE)
  
  # 宽度：每个面板 ~1.6 英寸
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
      # 轴线粗细
      axis.line.x = element_line(linewidth = 0.1),
      axis.line.y = element_line(linewidth = 0.1),
      # 刻度线粗细与长度
      axis.ticks.x = element_line(linewidth = 0.1),
      axis.ticks.y = element_line(linewidth = 0.1),
      axis.ticks.length = unit(4, "pt"),
      # 文字
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

## ===== 主流程：分别画 786_PCR 和 769_PCR =====
sheets <- tryCatch(readxl::excel_sheets(infile), error = function(e) character(0))

for (sn in c("786_PCR","769_PCR")){
  if (!sn %in% sheets){
    message("跳过【", sn, "】：sheet 不存在。")
    next
  }
  dat_long <- read_sheet_as_long(infile, sn)
  plot_facets_by_treatment(dat_long, sn)
}

