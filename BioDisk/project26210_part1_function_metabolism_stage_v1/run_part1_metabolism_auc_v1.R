#!/usr/bin/env Rscript

project <- "/mnt/BioDisk/ludingyu/project26210"
env_dir <- file.path(project, "envs", "part2_6_2_7_20260726_v1")
.libPaths(c(file.path(env_dir, "R_library"), .libPaths()))
suppressPackageStartupMessages({
  library(Matrix); library(data.table); library(AUCell); library(GSEABase); library(ggplot2)
})
set.seed(26210)

root <- file.path(project, "results", "part1_rebuild_v1")
input <- file.path(root, "07_pyscenic_v1", "input")
clone_stats <- file.path(root, "04_subclones_tymp", "cnv_clone_tymp_cnmf_association_v1", "patient_factor_high_clone_tymp_contrasts_v1.tsv")
out <- file.path(root, "09_metabolism_trajectory", "metabolism_v1")
invisible(lapply(file.path(out, c("figures", "tables", "source_data", "objects", "qc", "scFEA_input", "logs")), dir.create, recursive = TRUE, showWarnings = FALSE))

gate <- system2(file.path(project, "check_memory_gate.sh"), "part1_metabolism_auc_v1", stdout = file.path(out, "qc", "memory_gate_auc_v1.txt"), stderr = file.path(out, "qc", "memory_gate_auc_v1.txt"))
if (!identical(gate, 0L)) stop("Memory gate refused metabolism AUCell")

genes <- fread(file.path(input, "genes_v1.tsv"))$gene
cells <- fread(file.path(input, "cells_v1.tsv"))$exact_cell_id
meta <- fread(file.path(input, "cell_metadata_v1.tsv.gz"))
stopifnot(identical(meta$exact_cell_id, cells), length(cells) == 11950L)
con <- gzfile(file.path(input, "raw_counts_genes_by_cells_v1.mtx.gz"), "rb")
counts <- as(readMM(con), "dgCMatrix"); close(con)
rownames(counts) <- genes; colnames(counts) <- cells
stopifnot(identical(dim(counts), c(length(genes), length(cells))), length(counts@x) == 35713877L)

gmt <- file.path(env_dir, "sources", "scMetabolism", "data", "KEGG_metabolism_nc.gmt")
gene_sets <- GSEABase::getGmt(gmt)
targets <- c("Pyrimidine metabolism", "Cysteine and methionine metabolism", "Arginine and proline metabolism", "Glycerophospholipid metabolism")
stopifnot(all(targets %in% names(gene_sets)))
rankings <- AUCell::AUCell_buildRankings(counts, plotStats = FALSE, splitByBlocks = TRUE, BPPARAM = BiocParallel::SerialParam(progressbar = TRUE), verbose = TRUE)
auc <- AUCell::AUCell_calcAUC(gene_sets, rankings, aucMaxRank = ceiling(0.05 * nrow(counts)), nCores = 1L)
auc_matrix <- as.matrix(SummarizedExperiment::assay(auc))[, cells, drop = FALSE]
fwrite(as.data.table(t(auc_matrix), keep.rownames = "exact_cell_id"), file.path(out, "source_data", "KEGG_AUCell_all_pathways_v1.tsv.gz"), sep = "\t")

short <- c("Pyrimidine metabolism"="Pyrimidine", "Cysteine and methionine metabolism"="CysMet", "Arginine and proline metabolism"="ArgPro", "Glycerophospholipid metabolism"="Glycerophospholipid")
for (p in targets) meta[[paste0("AUC_", short[p])]] <- auc_matrix[p, ]
score_cols <- paste0("AUC_", unname(short))

contrasts <- fread(clone_stats)[distance_threshold == 0.25 & factor == "F1", .(curated_patient_id = infercnv_unit, factor_high_clone_id)]
stopifnot(nrow(contrasts) == 14L)
meta <- merge(meta, contrasts, by = "curated_patient_id", all.x = TRUE, sort = FALSE)
meta[, idx__ := match(exact_cell_id, cells)]; setorder(meta, idx__); meta[, idx__ := NULL]
meta[, clone_group := fifelse(clone_size_class_v1 == "primary_clone_ge20_v1" & formal_inferred_cnv_clone_v1 == factor_high_clone_id, "F1_high_primary_clone", fifelse(clone_size_class_v1 == "primary_clone_ge20_v1" & !is.na(factor_high_clone_id), "other_primary_clones", NA_character_))]
meta[, TYMP_group := fifelse(TYMP_raw_count_v1 > 0, "TYMP_positive", "TYMP_negative")]
saveRDS(meta, file.path(out, "objects", "metabolism_cell_metadata_four_auc_v1.rds"), compress = "xz")
fwrite(meta[, c("exact_cell_id", "study", "curated_patient_id", "TYMP_raw_count_v1", "TYMP_log1p_CP10K_v1", "F1", "clone_group", "TYMP_group", score_cols), with=FALSE], file.path(out, "source_data", "four_pathway_AUCell_cell_scores_v1.tsv.gz"), sep = "\t")

safe_wilcox <- function(x) if (length(x) < 2L || all(abs(x) < 1e-12)) if (length(x) >= 2L) 1 else NA_real_ else suppressWarnings(wilcox.test(x, mu=0, exact=FALSE)$p.value)
paired_stats <- function(group_col, low, high, contrast, min_cells=20L) {
  long <- melt(meta[!is.na(get(group_col))], id.vars=c("curated_patient_id", group_col), measure.vars=score_cols, variable.name="pathway", value.name="score")
  agg <- long[, .(median_score=median(score), n_cells=.N), by=c("curated_patient_id", "pathway", group_col)]
  setnames(agg, group_col, "group")
  wide <- dcast(agg[n_cells >= min_cells], curated_patient_id + pathway ~ group, value.var="median_score")
  if (!all(c(low,high) %in% names(wide))) return(data.table())
  wide <- wide[is.finite(get(low)) & is.finite(get(high))]
  detail <- wide[, .(curated_patient_id, pathway, low_score=get(low), high_score=get(high), delta=get(high)-get(low))]
  stat <- detail[, .(contrast=contrast, n_patients=.N, median_patient_delta=median(delta), fraction_patients_positive=mean(delta>0), paired_wilcox_PValue=safe_wilcox(delta)), by=pathway]
  stat[, FDR := p.adjust(paired_wilcox_PValue, "BH")]
  list(detail=detail[, contrast:=contrast], stat=stat)
}
tymp <- paired_stats("TYMP_group", "TYMP_negative", "TYMP_positive", "TYMP-positive_vs_negative", 20L)
clone <- paired_stats("clone_group", "other_primary_clones", "F1_high_primary_clone", "F1-high_clone_vs_other", 20L)
details <- rbindlist(list(tymp$detail, clone$detail), fill=TRUE); stats <- rbindlist(list(tymp$stat, clone$stat), fill=TRUE)
fwrite(details, file.path(out, "tables", "four_pathway_AUCell_patient_deltas_v1.tsv"), sep="\t")
fwrite(stats, file.path(out, "tables", "four_pathway_AUCell_paired_statistics_v1.tsv"), sep="\t")

cor_one <- function(predictor) {
  rbindlist(lapply(score_cols, function(sc) meta[, {ok <- is.finite(get(sc)) & is.finite(get(predictor)); .(n_cells=sum(ok), rho=if(sum(ok)>=30 && sd(get(sc)[ok])>0 && sd(get(predictor)[ok])>0) suppressWarnings(cor(get(sc)[ok],get(predictor)[ok],method="spearman")) else NA_real_)}, by=curated_patient_id][, `:=`(pathway=sc,predictor=predictor)]))
}
cors <- rbindlist(lapply(c("TYMP_log1p_CP10K_v1", "F1"), cor_one))[is.finite(rho)]
cor_summary <- cors[, .(n_patients=.N, median_patient_rho=median(rho), fraction_positive=mean(rho>0), signed_rank_PValue=safe_wilcox(rho)), by=.(pathway,predictor)]
cor_summary[, FDR := p.adjust(signed_rank_PValue,"BH"), by=predictor]
fwrite(cors, file.path(out,"tables","four_pathway_AUCell_correlations_by_patient_v1.tsv"),sep="\t")
fwrite(cor_summary, file.path(out,"tables","four_pathway_AUCell_correlation_summary_v1.tsv"),sep="\t")

stats[, pathway_label := sub("^AUC_", "", pathway)]; stats[, label := sprintf("%.3g%s", median_patient_delta, fifelse(FDR<0.05,"*",""))]
p1 <- ggplot(stats, aes(contrast, pathway_label, fill=median_patient_delta)) + geom_tile(colour="white") + geom_text(aes(label=label),size=3) + scale_fill_gradient2(low="#355F9D",mid="white",high="#C94755",midpoint=0) + theme_minimal(base_size=9) + theme(panel.grid=element_blank(),axis.text.x=element_text(angle=25,hjust=1)) + labs(x=NULL,y=NULL,fill="Median delta",title="Patient-paired KEGG/AUCell pathway differences")
ggsave(file.path(out,"figures","four_pathway_AUCell_paired_heatmap_v1.png"),p1,width=6.8,height=3.8,dpi=300,bg="white"); ggsave(file.path(out,"figures","four_pathway_AUCell_paired_heatmap_v1.pdf"),p1,width=6.8,height=3.8)
jsonlite::write_json(list(schema_version="part1_metabolism_auc_v1",completed_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),n_cells=length(cells),n_genes=length(genes),n_kegg_pathways=nrow(auc_matrix),n_target_pathways=4,method="KEGG gene sets with AUCell on raw UMI counts",statistical_unit="patient-paired medians",dependency_changes=FALSE), file.path(out,"qc","auc_metadata_v1.json"),auto_unbox=TRUE,pretty=TRUE)
writeLines("PASS", file.path(out,"AUC_COMPLETE")); cat("PART1_METABOLISM_AUC_V1_PASS\n")
