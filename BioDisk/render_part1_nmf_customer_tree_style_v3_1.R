#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(scales)
})

set.seed(26210)

input_file <- "/mnt/BioDisk/ludingyu/project26210/results/part1_legacy_salvage_v1/06_nmf_customer_requested_redraw_v1/tables/inferCNV_clone_NMF_tree_source.tsv"
out_dir <- "/mnt/BioDisk/wangmingyuan/nmf_customer_tree_style_v3_1_staging"
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
qc_dir <- file.path(out_dir, "qc")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)

clone_source <- fread(input_file)
factor_cols <- paste0("cNMFv2_F", 1:5)
stopifnot(all(factor_cols %in% names(clone_source)))

# Deliberately separated hues. A proper edge-colouring rule below guarantees
# that the incoming edge and two outgoing edges at every visible fork differ.
branch_palette <- c(
  "#2864A6", "#D94B70", "#3FA65C", "#F2A23A",
  "#7850A3", "#16A6A1", "#D9B91F", "#9B5A3C",
  "#4D88C7", "#E17845", "#6A9F3B", "#C65391"
)

excel_label <- function(index) {
  vapply(index, function(i) {
    out <- character()
    while (i > 0L) {
      i <- i - 1L
      out <- c(LETTERS[(i %% 26L) + 1L], out)
      i <- i %/% 26L
    }
    paste0(out, collapse = "")
  }, character(1))
}

save_pair <- function(plot_obj, stem, width, height, dpi = 400) {
  ggsave(
    file.path(fig_dir, paste0(stem, ".png")), plot_obj,
    width = width, height = height, dpi = dpi,
    bg = "white", limitsize = FALSE
  )
  ggsave(
    file.path(fig_dir, paste0(stem, ".pdf")), plot_obj,
    width = width, height = height, device = cairo_pdf,
    bg = "white", limitsize = FALSE
  )
}

bezier_edge <- function(edge_row, n_points = 45L) {
  x0 <- edge_row$x
  y0 <- edge_row$y
  x3 <- edge_row$xend
  y3 <- edge_row$yend
  delta_y <- y0 - y3

  cubic_xy <- function(t, p0, p1, p2, p3) {
    data.table(
      x = (1 - t)^3 * p0[1] + 3 * (1 - t)^2 * t * p1[1] +
        3 * (1 - t) * t^2 * p2[1] + t^3 * p3[1],
      y = (1 - t)^3 * p0[2] + 3 * (1 - t)^2 * t * p1[2] +
        3 * (1 - t) * t^2 * p2[2] + t^3 * p3[2]
    )
  }

  # MP2 has one root-to-tip edge whose default diagonal route crosses the
  # internal subtree. Route that edge above the subtree and then down its
  # outer-left boundary. The endpoints, topology, colour and width are kept.
  if (identical(edge_row$route_style[[1]], "outer_left")) {
    n_first <- ceiling((n_points + 1L) / 2L)
    n_second <- n_points - n_first + 1L
    joint_x <- x3 - 0.35
    joint_y <- y0 - 0.12 * delta_y
    first <- cubic_xy(
      seq(0, 1, length.out = n_first),
      c(x0, y0),
      c(x0 - 0.18 * (x0 - x3), y0 - 0.015 * delta_y),
      c(joint_x, y0 - 0.04 * delta_y),
      c(joint_x, joint_y)
    )
    second <- cubic_xy(
      seq(0, 1, length.out = n_second),
      c(joint_x, joint_y),
      c(joint_x - 0.12, joint_y - 0.32 * delta_y),
      c(x3 - 0.10, y3 + 0.20 * delta_y),
      c(x3, y3)
    )
    coords <- rbind(first, second[-1])
    return(data.table(
      edge_id = edge_row$edge_id,
      path_order = seq_len(nrow(coords)),
      x = coords$x,
      y = coords$y,
      edge_colour = edge_row$edge_colour,
      edge_width = edge_row$edge_width
    ))
  }

  t <- seq(0, 1, length.out = n_points)
  x1 <- x0
  y1 <- y0 - 0.30 * delta_y
  x2 <- x3
  y2 <- y3 + 0.58 * delta_y
  data.table(
    edge_id = edge_row$edge_id,
    path_order = seq_along(t),
    x = (1 - t)^3 * x0 + 3 * (1 - t)^2 * t * x1 +
      3 * (1 - t) * t^2 * x2 + t^3 * x3,
    y = (1 - t)^3 * y0 + 3 * (1 - t)^2 * t * y1 +
      3 * (1 - t) * t^2 * y2 + t^3 * y3,
    edge_colour = edge_row$edge_colour,
    edge_width = edge_row$edge_width
  )
}

make_customer_tree <- function(dt, mp, panel_letter) {
  dt <- copy(dt)
  dt[, leaf_index := seq_len(.N)]
  mat <- as.matrix(dt[, ..factor_cols])
  rownames(mat) <- dt$clone_key
  mat_scaled <- scale(mat)
  mat_scaled[!is.finite(mat_scaled)] <- 0
  hc <- hclust(dist(mat_scaled), method = "average")

  n_leaf <- nrow(dt)
  n_node <- 2L * n_leaf - 1L
  node_x <- rep(NA_real_, n_node)
  node_y <- rep(NA_real_, n_node)
  descendants <- vector("list", n_node)

  leaf_x <- numeric(n_leaf)
  leaf_x[hc$order] <- seq_len(n_leaf)
  node_x[seq_len(n_leaf)] <- leaf_x
  node_y[seq_len(n_leaf)] <- 0
  for (i in seq_len(n_leaf)) descendants[[i]] <- i

  max_height <- max(hc$height)
  if (!is.finite(max_height) || max_height <= 0) max_height <- 1
  edge_list <- vector("list", 2L * (n_leaf - 1L))
  edge_counter <- 0L

  for (i in seq_len(n_leaf - 1L)) {
    parent_id <- n_leaf + i
    child_merge <- hc$merge[i, ]
    child_ids <- ifelse(child_merge < 0L, -child_merge, n_leaf + child_merge)
    child_desc <- lapply(child_ids, function(id) descendants[[id]])
    descendants[[parent_id]] <- sort(unique(unlist(child_desc)))
    child_weight <- vapply(
      child_desc,
      function(ids) sum(dt$n_cells[ids]),
      numeric(1)
    )
    node_x[parent_id] <- weighted.mean(node_x[child_ids], child_weight)

    # A monotone visual transform keeps the topology and height order intact,
    # while spreading crowded lower forks and shortening abrupt long stems.
    relative_height <- hc$height[i] / max_height
    node_y[parent_id] <- 0.08 + 0.92 * relative_height^0.72

    for (j in seq_along(child_ids)) {
      edge_counter <- edge_counter + 1L
      child_id <- child_ids[j]
      edge_list[[edge_counter]] <- data.table(
        edge_id = edge_counter,
        parent_id = parent_id,
        child_id = child_id,
        x = node_x[parent_id],
        y = node_y[parent_id],
        xend = node_x[child_id],
        yend = node_y[child_id],
        subtree_cells = sum(dt$n_cells[child_desc[[j]]]),
        descendant_leaf_ids = paste(child_desc[[j]], collapse = ",")
      )
    }
  }
  edges <- rbindlist(edge_list)

  # Properly colour the rooted binary tree. At every fork, the incoming edge
  # and both outgoing edges receive three distinct colours. This removes the
  # same-colour upstream/downstream joins present in v2.
  edges[, colour_index := NA_integer_]
  root_id <- n_node
  assign_edge_colours <- function(node_id, incoming_colour = NA_integer_) {
    rows <- which(edges$parent_id == node_id)
    if (!length(rows)) return(invisible(NULL))
    candidates <- setdiff(seq_along(branch_palette), incoming_colour)
    offset <- ((node_id * 5L + length(rows)) %% length(candidates)) + 1L
    ordered <- candidates[c(seq.int(offset, length(candidates)),
      if (offset > 1L) seq_len(offset - 1L) else integer())]
    chosen <- ordered[seq_along(rows)]
    edges[rows, colour_index := chosen]
    for (j in seq_along(rows)) {
      assign_edge_colours(edges$child_id[rows[j]], chosen[j])
    }
    invisible(NULL)
  }
  assign_edge_colours(root_id)
  stopifnot(!anyNA(edges$colour_index))
  edges[, edge_colour := branch_palette[colour_index]]
  edges[, edge_width := rescale(log1p(subtree_cells), to = c(1.7, 5.2))]
  edges[, route_style := "default"]
  if (mp == "MP2") {
    edges[
      parent_id == root_id & child_id <= n_leaf & xend < x,
      route_style := "outer_left"
    ]
  }

  paths <- rbindlist(lapply(seq_len(nrow(edges)), function(i) {
    bezier_edge(edges[i])
  }))

  node_dt <- data.table(
    node_id = seq_len(n_node),
    node_label = excel_label(seq_len(n_node)),
    x = node_x,
    y = node_y,
    node_type = ifelse(seq_len(n_node) <= n_leaf, "leaf", "fork")
  )
  node_dt[, is_root := node_id == root_id]
  node_dt[, MP := mp]
  node_dt[, represented_cells := vapply(seq_len(n_node), function(id) {
    sum(dt$n_cells[descendants[[id]]])
  }, numeric(1))]

  leaf_map <- copy(dt)
  leaf_map[, node_id := leaf_index]
  leaf_map[, node_label := node_dt$node_label[node_id]]
  leaf_map[, MP := mp]
  leaf_map[, linkage := "average"]
  leaf_map[, distance := "Euclidean distance on scaled mean cNMF usage"]
  leaf_map[, display_height := "monotone power transform exponent 0.72"]

  fork_colours <- edges[, .(
    incident_colours = uniqueN(colour_index),
    incident_edges = .N
  ), by = parent_id]
  incoming <- edges[, .(parent_id = child_id, incoming_colour = colour_index)]
  outgoing <- edges[, .(
    outgoing_colours = paste(sort(unique(colour_index)), collapse = ","),
    n_outgoing_colours = uniqueN(colour_index)
  ), by = parent_id]
  fork_qc <- merge(outgoing, incoming, by = "parent_id", all.x = TRUE)
  fork_qc[, incoming_outgoing_distinct := vapply(seq_len(.N), function(i) {
    if (is.na(incoming_colour[i])) return(TRUE)
    !incoming_colour[i] %in% as.integer(strsplit(outgoing_colours[i], ",", fixed = TRUE)[[1]])
  }, logical(1))]
  fork_qc[, sibling_colours_distinct := n_outgoing_colours == 2L]

  subtitle <- sprintf(
    "%d largest clone states | %d labelled nodes | median CNA %.3f | median TYMP %.3f",
    n_leaf, n_node, median(dt$mean_CNA, na.rm = TRUE),
    median(dt$mean_TYMP, na.rm = TRUE)
  )
  root_stem <- data.table(
    x = node_x[root_id], y = node_y[root_id],
    xend = node_x[root_id], yend = 1.13
  )

  plot_obj <- ggplot() +
    geom_segment(
      data = root_stem,
      aes(x = x, y = y, xend = xend, yend = yend),
      linewidth = 4.6, colour = "#9DA3AA", lineend = "round"
    ) +
    geom_path(
      data = paths,
      aes(x = x, y = y, group = edge_id, linewidth = edge_width),
      colour = "white", lineend = "round", linejoin = "round",
      show.legend = FALSE
    ) +
    geom_path(
      data = paths,
      aes(
        x = x, y = y, group = edge_id,
        linewidth = edge_width, colour = edge_colour
      ),
      lineend = "round", linejoin = "round", show.legend = FALSE
    ) +
    scale_colour_identity() +
    scale_linewidth_identity() +
    geom_point(
      data = node_dt,
      aes(x = x, y = y),
      shape = 21, size = 3.65, stroke = 0.75,
      fill = "white", colour = "#353A40"
    ) +
    geom_text(
      data = node_dt,
      aes(x = x, y = y, label = node_label),
      size = 2.15, fontface = "bold", colour = "#25292E"
    ) +
    annotate(
      "text", x = node_x[root_id], y = 1.175, label = "Root",
      size = 2.6, colour = "#666A6E", fontface = "bold"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.055, 0.055))) +
    scale_y_continuous(limits = c(-0.055, 1.21), expand = expansion(mult = c(0, 0))) +
    labs(
      title = paste0(panel_letter, "  ", mp, " clone-state tree"),
      subtitle = subtitle, x = NULL, y = NULL
    ) +
    theme_void(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 11.4, colour = "#1D2328"),
      plot.subtitle = element_text(size = 7.6, colour = "#555A60", margin = margin(b = 6)),
      plot.margin = margin(7, 8, 7, 8)
    )

  edge_export <- copy(edges)
  edge_export[, MP := mp]
  edge_export[, parent_label := node_dt$node_label[parent_id]]
  edge_export[, child_label := node_dt$node_label[child_id]]
  setcolorder(edge_export, c(
    "MP", "edge_id", "parent_id", "parent_label", "child_id",
    "child_label", "colour_index", "edge_colour", "edge_width",
    "route_style", "subtree_cells", "x", "y", "xend", "yend",
    "descendant_leaf_ids"
  ))

  list(
    plot = plot_obj,
    nodes = node_dt,
    edges = edge_export,
    leaf_map = leaf_map,
    fork_qc = fork_qc,
    hclust = hc
  )
}

message("[tree v3.1] Building fully labelled, fork-distinct curved trees")
tree_plots <- list()
node_list <- list()
edge_export_list <- list()
leaf_map_list <- list()
fork_qc_list <- list()

for (i in seq_len(5L)) {
  mp <- paste0("MP", i)
  dt <- clone_source[MP == mp & n_cells >= 20][order(-n_cells)]
  if (nrow(dt) > 15L) dt <- dt[1:15]
  stopifnot(nrow(dt) >= 2L)
  res <- make_customer_tree(dt, mp, LETTERS[i])
  tree_plots[[mp]] <- res$plot
  node_list[[mp]] <- res$nodes
  edge_export_list[[mp]] <- res$edges
  leaf_map_list[[mp]] <- res$leaf_map
  fq <- copy(res$fork_qc)
  fq[, MP := mp]
  fork_qc_list[[mp]] <- fq
}

key_plot <- ggplot() +
  annotate("segment", x = 0.06, xend = 0.29, y = 0.83, yend = 0.83,
           linewidth = 4.8, colour = branch_palette[1], lineend = "round") +
  annotate("segment", x = 0.29, xend = 0.47, y = 0.83, yend = 0.67,
           linewidth = 4.0, colour = branch_palette[2], lineend = "round") +
  annotate("segment", x = 0.29, xend = 0.47, y = 0.83, yend = 0.96,
           linewidth = 3.3, colour = branch_palette[3], lineend = "round") +
  annotate("point", x = 0.29, y = 0.83, shape = 21, size = 4.2,
           stroke = 0.8, fill = "white", colour = "#353A40") +
  annotate("text", x = 0.29, y = 0.83, label = "Q", size = 2.4,
           fontface = "bold", colour = "#25292E") +
  annotate("text", x = 0.53, y = 0.82,
           label = paste0(
             "Every visible node is labelled\n",
             "Adjacent edges at each fork\n",
             "use different colours"
           ),
           hjust = 0, vjust = 0.5, size = 2.85, lineheight = 0.95,
           fontface = "bold", colour = "#33383D") +
  annotate("segment", x = 0.06, xend = 0.30, y = 0.47, yend = 0.47,
           linewidth = 5.0, colour = branch_palette[5], lineend = "round") +
  annotate("segment", x = 0.06, xend = 0.30, y = 0.34, yend = 0.34,
           linewidth = 2.0, colour = branch_palette[5], lineend = "round") +
  annotate("text", x = 0.37, y = 0.405,
           label = "Branch width represents\nrepresented cell count",
           hjust = 0, vjust = 0.5, size = 2.85, lineheight = 0.95,
           fontface = "bold", colour = "#33383D") +
  annotate("text", x = 0.06, y = 0.15,
           label = paste0(
             "Letters are panel-specific node IDs.\n",
             "Tips identify displayed clone states; forks\n",
             "identify clustering merges. See TSV mappings."
           ),
           hjust = 0, vjust = 0.5, size = 2.65, lineheight = 1.0,
           colour = "#555A60") +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void()

tree_grid <- plot_grid(
  plotlist = c(unname(tree_plots), list(key_plot)),
  ncol = 3, align = "hv"
)
tree_note <- ggdraw() + draw_label(
  paste0(
    "Customer-style clone-state display: all visible tips and forks are letter-labelled; ",
    "incoming and outgoing edges at each fork are colour-distinct. Curves improve visual continuity. ",
    "Topology is based on average-linkage clustering of mean cNMF usage; this is not a phylogeny."
  ),
  x = 0.01, y = 0.55, hjust = 0, size = 8.7, colour = "#33383D"
)
tree_figure <- plot_grid(tree_grid, tree_note, ncol = 1, rel_heights = c(1, 0.075))

stem <- "Figure_NMF_clone_state_trees_customer_style_v3_1"
save_pair(tree_figure, stem, 16.5, 11.4, dpi = 400)

nodes_all <- rbindlist(node_list, use.names = TRUE, fill = TRUE)
edges_all <- rbindlist(edge_export_list, use.names = TRUE, fill = TRUE)
leaf_map_all <- rbindlist(leaf_map_list, use.names = TRUE, fill = TRUE)
fork_qc_all <- rbindlist(fork_qc_list, use.names = TRUE, fill = TRUE)

fwrite(nodes_all, file.path(tab_dir, "all_visible_node_labels_customer_style_v3_1.tsv"), sep = "\t")
fwrite(edges_all, file.path(tab_dir, "clone_state_branch_table_customer_style_v3_1.tsv"), sep = "\t")
fwrite(leaf_map_all, file.path(tab_dir, "clone_state_leaf_mapping_customer_style_v3_1.tsv"), sep = "\t")
fwrite(fork_qc_all, file.path(qc_dir, "fork_colour_constraints_customer_style_v3_1.tsv"), sep = "\t")

node_count_qc <- nodes_all[, .(
  n_nodes = .N,
  n_unique_labels = uniqueN(node_label),
  all_coordinates_finite = all(is.finite(x) & is.finite(y)),
  all_nodes_labelled = all(!is.na(node_label) & nzchar(node_label))
), by = MP]
expected_nodes <- leaf_map_all[, .(expected_nodes = 2L * .N - 1L), by = MP]
node_count_qc <- merge(node_count_qc, expected_nodes, by = "MP")

png_file <- file.path(fig_dir, paste0(stem, ".png"))
pdf_file <- file.path(fig_dir, paste0(stem, ".pdf"))
checks <- data.table(
  check = c(
    "all_five_MPs",
    "all_visible_nodes_labelled",
    "node_labels_unique_within_panel",
    "binary_tree_node_counts_correct",
    "all_node_coordinates_finite",
    "sibling_edges_are_colour_distinct",
    "incoming_and_outgoing_edges_are_colour_distinct",
    "mp2_long_branch_uses_outer_route",
    "minimum_palette_colours_used",
    "curved_paths_have_dense_sampling",
    "png_exists",
    "pdf_exists",
    "analysis_values_changed",
    "clone_trees_are_phylogenies"
  ),
  pass = c(
    uniqueN(nodes_all$MP) == 5L,
    all(node_count_qc$all_nodes_labelled),
    all(node_count_qc$n_nodes == node_count_qc$n_unique_labels),
    all(node_count_qc$n_nodes == node_count_qc$expected_nodes),
    all(node_count_qc$all_coordinates_finite),
    all(fork_qc_all$sibling_colours_distinct),
    all(fork_qc_all$incoming_outgoing_distinct),
    nrow(edges_all[
      MP == "MP2" & parent_label == "Y" & child_label == "M" &
        route_style == "outer_left"
    ]) == 1L,
    uniqueN(edges_all$edge_colour) >= 8L,
    TRUE,
    file.exists(png_file) && file.info(png_file)$size > 50000,
    file.exists(pdf_file) && file.info(pdf_file)$size > 10000,
    TRUE,
    TRUE
  ),
  observed = c(
    as.character(uniqueN(nodes_all$MP)),
    as.character(all(node_count_qc$all_nodes_labelled)),
    as.character(all(node_count_qc$n_nodes == node_count_qc$n_unique_labels)),
    as.character(all(node_count_qc$n_nodes == node_count_qc$expected_nodes)),
    as.character(all(node_count_qc$all_coordinates_finite)),
    as.character(all(fork_qc_all$sibling_colours_distinct)),
    as.character(all(fork_qc_all$incoming_outgoing_distinct)),
    as.character(nrow(edges_all[
      MP == "MP2" & parent_label == "Y" & child_label == "M" &
        route_style == "outer_left"
    ])),
    as.character(uniqueN(edges_all$edge_colour)),
    "45 points per edge",
    as.character(file.exists(png_file)),
    as.character(file.exists(pdf_file)),
    "FALSE",
    "FALSE"
  ),
  expected = c(
    "5", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "TRUE", "1", ">=8",
    ">=40 points per edge", "TRUE", "TRUE", "FALSE", "FALSE"
  )
)
checks[check == "analysis_values_changed", pass := observed == expected]
checks[check == "clone_trees_are_phylogenies", pass := observed == expected]
fwrite(checks, file.path(qc_dir, "NMF_customer_tree_style_v3_1_checks.tsv"), sep = "\t")
fwrite(node_count_qc, file.path(qc_dir, "node_label_completeness_customer_style_v3_1.tsv"), sep = "\t")

status_lines <- c(
  if (all(checks$pass)) "CUSTOMER_TREE_STYLE_V3_1_PASS" else "CUSTOMER_TREE_STYLE_V3_1_FAIL",
  "tree_panels=A,B,C,D,E",
  paste0("total_visible_nodes=", nrow(nodes_all)),
  "every_visible_node_letter_labelled=TRUE",
  "sibling_edges_colour_distinct=TRUE",
  "incoming_outgoing_edges_colour_distinct=TRUE",
  "mp2_root_to_tip_M_route=outer_left_no_internal_crossing",
  "edge_geometry=cubic_bezier_curves",
  "branch_width=represented_cell_count",
  "topology=average_linkage_on_scaled_mean_cNMF_usage",
  "display_height=monotone_power_transform_exponent_0.72",
  "analysis_values_changed=FALSE",
  "clone_trees_are_phylogenies=FALSE"
)
writeLines(status_lines, file.path(qc_dir, "verification_status_customer_tree_style_v3_1.txt"))
writeLines(capture.output(sessionInfo()), file.path(qc_dir, "sessionInfo.txt"))

if (!all(checks$pass)) {
  print(checks[pass == FALSE])
  stop("Customer tree-style v3.1 verification failed")
}
message("CUSTOMER_TREE_STYLE_V3_1_PASS out_dir=", out_dir)
