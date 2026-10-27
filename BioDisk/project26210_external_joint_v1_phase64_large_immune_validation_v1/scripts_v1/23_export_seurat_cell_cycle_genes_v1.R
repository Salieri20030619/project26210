args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: Rscript 23_export_seurat_cell_cycle_genes_v1.R <output.tsv>")

suppressPackageStartupMessages(library(Seurat))
data("cc.genes.updated.2019", package = "Seurat")

genes <- data.frame(
  phase = c(rep("S", length(cc.genes.updated.2019$s.genes)),
            rep("G2M", length(cc.genes.updated.2019$g2m.genes))),
  gene = c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes),
  source = "Seurat::cc.genes.updated.2019",
  stringsAsFactors = FALSE
)
write.table(genes, args[[1L]], sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Seurat=%s\nS_genes=%d\nG2M_genes=%d\n",
            as.character(packageVersion("Seurat")),
            sum(genes$phase == "S"), sum(genes$phase == "G2M")))
