args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: Rscript script ligand_target_matrix.rds lr_network.rds output_dir")
ligand_target_path <- args[[1]]
lr_network_path <- args[[2]]
out <- args[[3]]
if (dir.exists(out)) stop("Refusing to overwrite output directory")
dir.create(out, recursive = TRUE)

suppressPackageStartupMessages(library(CellChat))
suppressPackageStartupMessages(library(nichenetr))
data("CellChatDB.human", package = "CellChat")
ligand_target_matrix <- readRDS(ligand_target_path)
lr_network <- readRDS(lr_network_path)

genes <- c("MST1", "MST1R", "ICAM5", "ITGAM", "ITGB2", "TYMP", "HIF1A")
lt <- data.frame(
  gene = genes,
  in_ligand_target_rows = genes %in% rownames(ligand_target_matrix),
  in_ligand_target_columns = genes %in% colnames(ligand_target_matrix),
  stringsAsFactors = FALSE
)

lr <- as.data.frame(lr_network)
names(lr)[1:2] <- c("from", "to")
targeted_lr <- lr[lr$from %in% c("MST1", "ICAM5") | lr$to %in% c("MST1R", "ITGAM", "ITGB2"), , drop = FALSE]

cc <- CellChatDB.human$interaction
cc_cols <- intersect(c("interaction_name", "interaction_name_2", "pathway_name", "ligand", "receptor"), names(cc))
targeted_cc <- cc[apply(cc, 1, function(z) any(grepl("MST1|MST1R|ICAM5|ITGAM|ITGB2", as.character(z)))), cc_cols, drop = FALSE]

write.table(lt, file.path(out, "nichenet_gene_coverage_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(targeted_lr, file.path(out, "nichenet_targeted_lr_network_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(targeted_cc, file.path(out, "cellchat_targeted_interactions_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
status <- data.frame(
  metric = c("cellchat_interactions", "cellchat_targeted_rows", "nichenet_lr_rows", "nichenet_targeted_lr_rows", "ligand_target_rows", "ligand_target_columns"),
  value = c(nrow(cc), nrow(targeted_cc), nrow(lr), nrow(targeted_lr), nrow(ligand_target_matrix), ncol(ligand_target_matrix))
)
write.table(status, file.path(out, "resource_dimensions_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
