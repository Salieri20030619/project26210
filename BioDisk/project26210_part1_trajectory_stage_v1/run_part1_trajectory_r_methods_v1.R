#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(monocle3)
  library(SingleCellExperiment)
  library(slingshot)
  library(ggplot2)
})

set.seed(26210)
options(mc.cores = 8L)
project <- "/mnt/BioDisk/ludingyu/project26210"
out <- file.path(project, "results/part1_rebuild_v1/09_metabolism_trajectory/trajectory_v1")
inp <- file.path(out, "input")
qc <- file.path(out, "qc")
dir.create(file.path(out, "figures"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "tables"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "source_data"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(out, "objects"), recursive=TRUE, showWarnings=FALSE)
stopifnot(file.exists(file.path(inp, "INPUT_COMPLETE")), !file.exists(file.path(out, "R_METHODS_COMPLETE")))

memory_gate <- function(stage) {
  ans <- system2(file.path(project, "check_memory_gate.sh"), stage, stdout=TRUE, stderr=TRUE)
  status <- attr(ans, "status"); if (is.null(status)) status <- 0L
  cat(stage, status, paste(ans, collapse=" | "), sep="\t", file=file.path(qc, "memory_gate_log_v1.tsv"), append=TRUE)
  cat("\n", file=file.path(qc, "memory_gate_log_v1.tsv"), append=TRUE)
  if (status != 0L) stop("Memory gate blocked ", stage)
}
scale01 <- function(x) {
  x <- as.numeric(x); ok <- is.finite(x); out <- rep(NA_real_, length(x))
  if (!any(ok)) return(out); r <- range(x[ok]); out[ok] <- if (diff(r)==0) 0.5 else (x[ok]-r[1])/diff(r); out
}
safe_cor <- function(x,y) {
  ok <- is.finite(x)&is.finite(y)
  if(sum(ok)<20L || sd(x[ok])==0 || sd(y[ok])==0) return(NA_real_)
  suppressWarnings(cor(x[ok],y[ok],method="spearman"))
}
safe_wilcox <- function(x) {
  x <- x[is.finite(x)]; if(length(x)<2L) return(NA_real_)
  if(all(abs(x)<1e-12)) return(1)
  suppressWarnings(wilcox.test(x,mu=0,exact=FALSE)$p.value)
}

memory_gate("part1_trajectory_r_input_v1")
cells <- fread(file.path(inp,"cells_v1.tsv"),header=FALSE)$V1
genes <- fread(file.path(inp,"genes_v1.tsv"),header=FALSE)$V1
md <- fread(file.path(inp,"cell_metadata_v1.tsv.gz")); stopifnot(identical(md$exact_cell_id,cells), nrow(md)==11950L)
counts <- readMM(gzfile(file.path(inp,"raw_counts_genes_by_cells_v1.mtx.gz")))
counts <- as(counts,"dgCMatrix"); rownames(counts)<-genes; colnames(counts)<-cells
stopifnot(nrow(counts)==2500L,ncol(counts)==11950L,all(counts@x>=0),all(abs(counts@x-round(counts@x))<1e-8))
lat_dt <- fread(file.path(inp,"scvi_latent_v1.tsv.gz"),header=TRUE,check.names=FALSE); umap_dt <- fread(file.path(inp,"scvi_umap_v1.tsv.gz"),header=TRUE,check.names=FALSE)
stopifnot(identical(lat_dt[[1]],cells),identical(umap_dt[[1]],cells))
lat <- as.matrix(lat_dt[,-1]); umap <- as.matrix(umap_dt[,-1]); rownames(lat)<-cells; rownames(umap)<-cells
colnames(umap)<-c("UMAP1","UMAP2")

cell_md <- as.data.frame(md); rownames(cell_md)<-cell_md$exact_cell_id
gene_md <- data.frame(gene_short_name=genes,row.names=genes)
cds <- new_cell_data_set(counts,cell_metadata=cell_md,gene_metadata=gene_md)
SingleCellExperiment::reducedDims(cds)$PCA <- lat[,seq_len(min(30L,ncol(lat))),drop=FALSE]
SingleCellExperiment::reducedDims(cds)$UMAP <- umap

memory_gate("part1_trajectory_monocle3_graph_v1")
cds <- cluster_cells(cds,reduction_method="UMAP",k=30L,cluster_method="leiden",num_iter=2L,random_seed=26210L,verbose=TRUE)
cds <- learn_graph(cds,use_partition=TRUE,close_loop=TRUE,learn_graph_control=list(minimal_branch_len=10L,euclidean_distance_ratio=1),verbose=TRUE)
part <- as.character(partitions(cds,reduction_method="UMAP")[cells])
root_cells <- md[root_candidate_v1==TRUE,exact_cell_id]
partition_roots <- md[,.(exact_cell_id=exact_cell_id[which.min(progression_anchor_v1)]),by=.(partition=part)]$exact_cell_id
root_cells <- unique(c(root_cells,partition_roots))
cds <- order_cells(cds,reduction_method="UMAP",root_cells=root_cells)
pt <- pseudotime(cds,reduction_method="UMAP")[cells]
md[,Monocle3_pseudotime_v1:=scale01(pt)]
md[,Monocle3_pseudotime_patient_v1:=scale01(Monocle3_pseudotime_v1),by=curated_patient_id]
md[,`:=`(Monocle3_cluster_v1=as.character(clusters(cds,reduction_method="UMAP")[cells]),Monocle3_partition_v1=part)]
if(sum(is.finite(md$Monocle3_pseudotime_v1)) < 0.95*nrow(md)) stop("Too few finite Monocle3 pseudotimes")

pg <- principal_graph(cds)[["UMAP"]]; pgc <- t(principal_graph_aux(cds)[["UMAP"]]$dp_mst)
pge <- igraph::as_data_frame(pg,what="edges")
segments <- data.table(from=pge$from,to=pge$to,x=pgc[pge$from,1],y=pgc[pge$from,2],xend=pgc[pge$to,1],yend=pgc[pge$to,2])
fwrite(segments,file.path(out,"source_data/monocle3_principal_graph_segments_v1.tsv"),sep="\t")

memory_gate("part1_trajectory_monocle3_graph_test_v1")
gt <- as.data.table(graph_test(cds,neighbor_graph="principal_graph",cores=8L),keep.rownames="gene")
setorder(gt,q_value)
fwrite(gt,file.path(out,"tables/monocle3_graph_test_hvg_v1.tsv.gz"),sep="\t")

memory_gate("part1_trajectory_slingshot_v1")
sce <- SingleCellExperiment(assays=list(counts=counts),colData=S4Vectors::DataFrame(cell_md[cells,,drop=FALSE]))
reducedDims(sce)$PCA <- lat[,seq_len(min(30L,ncol(lat))),drop=FALSE]
colData(sce)$trajectory_cluster <- factor(md$Monocle3_cluster_v1)
root_cluster <- md[,.(anchor=median(progression_anchor_v1)),by=Monocle3_cluster_v1][which.min(anchor),Monocle3_cluster_v1]
sce <- slingshot(sce,clusterLabels="trajectory_cluster",reducedDim="PCA",start.clus=root_cluster,shrink=TRUE,extend="y",reweight=TRUE,reassign=TRUE)
sl <- slingPseudotime(sce,na=FALSE); if(is.null(dim(sl))) sl<-matrix(sl,ncol=1L)
for(j in seq_len(ncol(sl))){rho_j<-safe_cor(sl[,j],md$progression_anchor_v1);if(is.finite(rho_j)&&rho_j<0) sl[,j]<-max(sl[,j])-sl[,j];sl[,j]<-scale01(sl[,j])}
lineage_qc <- data.table(lineage=seq_len(ncol(sl)),n_finite=colSums(is.finite(sl)),rho_anchor=vapply(seq_len(ncol(sl)),function(j)safe_cor(sl[,j],md$progression_anchor_v1),numeric(1)))
best <- lineage_qc[order(-n_finite,-rho_anchor),lineage][1]
md[,Slingshot_pseudotime_v1:=sl[,best]]
md[,Slingshot_pseudotime_patient_v1:=scale01(Slingshot_pseudotime_v1),by=curated_patient_id]
fwrite(data.table(exact_cell_id=cells,sl),file.path(out,"source_data/slingshot_all_lineages_v1.tsv.gz"),sep="\t")
fwrite(lineage_qc[,`:=`(selected=lineage==best,root_cluster=root_cluster)],file.path(out,"tables/slingshot_lineage_qc_v1.tsv"),sep="\t")

targets <- c("progression_anchor_v1","TYMP_log1p_CP10K_v1","F1","infercnv_burden_over_threshold_v1","AUC_Pyrimidine","AUC_CysMet","AUC_ArgPro","AUC_Glycerophospholipid","Pyrimidine_flux","CysMet_flux","ArgPro_flux","Choline_proxy")
methods <- c("Monocle3_pseudotime_patient_v1","Slingshot_pseudotime_patient_v1")
corr <- rbindlist(lapply(methods,function(method)rbindlist(lapply(targets,function(target)md[,.(n_cells=.N,rho=safe_cor(get(method),get(target))),by=curated_patient_id][,`:=`(method=method,target=target)]))))
corr <- corr[is.finite(rho)]
summary <- corr[,.(n_patients=.N,median_patient_rho=median(rho),fraction_positive=mean(rho>0),signed_rank_PValue=safe_wilcox(rho)),by=.(method,target)]
summary[,FDR:=p.adjust(signed_rank_PValue,"BH"),by=method]
fwrite(corr,file.path(out,"tables/trajectory_r_patient_correlations_v1.tsv"),sep="\t")
fwrite(summary,file.path(out,"tables/trajectory_r_correlation_summary_v1.tsv"),sep="\t")
fwrite(md,file.path(out,"source_data/trajectory_r_cell_scores_v1.tsv.gz"),sep="\t")

plot_dt <- md[,.(UMAP1=umap[,1],UMAP2=umap[,2],Monocle3=Monocle3_pseudotime_v1,Slingshot=Slingshot_pseudotime_v1)]
plot_one <- function(value,title,stem){
  p<-ggplot(plot_dt,aes(UMAP1,UMAP2,color=.data[[value]]))+geom_point(size=.18,alpha=.8)+scale_color_viridis_c(option="plasma",na.value="#D9D9D9")+coord_equal()+theme_void(base_size=9)+labs(title=title,color="Pseudotime")
  ggsave(file.path(out,"figures",paste0(stem,".png")),p,width=4.4,height=3.8,dpi=400,bg="white")
  ggsave(file.path(out,"figures",paste0(stem,".pdf")),p,width=4.4,height=3.8)
}
plot_one("Monocle3","Monocle3 cross-sectional pseudotime","Monocle3_pseudotime_v1")
plot_one("Slingshot","Slingshot sensitivity pseudotime","Slingshot_pseudotime_v1")
saveRDS(cds,file.path(out,"objects/monocle3_cds_v1.rds"),compress="xz")
writeLines(capture.output(sessionInfo()),file.path(out,"qc/R_session_info_v1.txt"))
metadata <- list(schema_version="part1_trajectory_r_methods_v1",completed_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),n_cells=nrow(md),n_genes=nrow(counts),n_monocle_partitions=length(unique(part)),n_monocle_clusters=length(unique(md$Monocle3_cluster_v1)),n_finite_monocle=sum(is.finite(md$Monocle3_pseudotime_v1)),n_slingshot_lineages=ncol(sl),selected_slingshot_lineage=best,n_graph_test_fdr_lt_0_05=sum(gt$q_value<0.05,na.rm=TRUE),root_anchor="independent of TYMP,F1,metabolism",interpretation="cross-sectional state ordering;not lineage tracing or longitudinal time")
jsonlite::write_json(metadata,file.path(out,"qc/trajectory_r_metadata_v1.json"),auto_unbox=TRUE,pretty=TRUE)
writeLines(format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),file.path(out,"R_METHODS_COMPLETE"))
cat("PART1_TRAJECTORY_R_METHODS_V1_PASS\n")
