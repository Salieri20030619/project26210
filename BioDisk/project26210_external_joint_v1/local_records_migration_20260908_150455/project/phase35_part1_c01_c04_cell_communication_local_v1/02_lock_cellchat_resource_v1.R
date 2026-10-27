args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript script output_dir")
out <- args[[1]]
if (dir.exists(out)) stop("Refusing to overwrite output directory")
dir.create(out, recursive = TRUE)
suppressPackageStartupMessages(library(CellChat))
data("CellChatDB.human", package = "CellChat")

saveRDS(CellChatDB.human, file.path(out, "CellChatDB.human_v1.rds"), version = 3)
interaction <- CellChatDB.human$interaction
complex <- CellChatDB.human$complex
cofactor <- CellChatDB.human$cofactor
simple <- unique(c(as.character(interaction$ligand), as.character(interaction$receptor)))
complex_members <- unique(as.character(unlist(complex, use.names = FALSE)))
cofactor_members <- unique(as.character(unlist(cofactor, use.names = FALSE)))
genes <- sort(unique(c(simple, complex_members, cofactor_members)))
genes <- genes[!is.na(genes) & nzchar(genes)]
write.table(data.frame(gene = genes), file.path(out, "cellchat_required_genes_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(interaction, file.path(out, "cellchat_interactions_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
status <- data.frame(
  cellchat_version = as.character(packageVersion("CellChat")),
  n_interactions = nrow(interaction),
  n_complexes = nrow(complex),
  n_cofactors = nrow(cofactor),
  n_exported_gene_tokens = length(genes)
)
write.table(status, file.path(out, "status_v1.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
