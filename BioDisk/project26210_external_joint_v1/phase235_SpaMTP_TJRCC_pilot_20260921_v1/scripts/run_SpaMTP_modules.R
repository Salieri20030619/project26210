options(width=160,Seurat.object.assay.version="v5")
suppressPackageStartupMessages({library(Seurat);library(SeuratObject);library(Matrix);library(data.table);library(jsonlite)})
P <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase235_SpaMTP_TJRCC_pilot_20260921_v1"
for (name in c("helperFunctions.R","LoadData.R","SMPreprocessing.R","DimensionalityReduction.R")) source(file.path(P,"metadata",name),local=.GlobalEnv)
stopifnot(!requireNamespace("SpaMTP",quietly=TRUE))
dir.create(file.path(P,"results","SpaMTP_objects"),showWarnings=FALSE)
allqc <- list()
run_one <- function(label, binned=FALSE) {
  suffix <- if(binned) "_bin4_SpaMTP_input.csv" else "_SpaMTP_input.csv"
  path <- file.path(P,"results","matrices",paste0(label,suffix))
  cat("IMPORT",label,binned,"\n");flush.console()
  obj <- ReadSM_mtx(path,assay="SPM",project.name=label,verbose=FALSE)
  input <- as.data.frame(data.table::fread(path,check.names=FALSE))
  expected <- t(as.matrix(input[,setdiff(colnames(input),c("x","y")),drop=FALSE]))
  counts <- as.matrix(LayerData(obj,assay="SPM",layer="counts"))
  stopifnot(all(dim(counts)==dim(expected)),max(abs(counts-expected))==0)
  xy <- GetTissueCoordinates(obj,image="fov")
  stopifnot(identical(as.character(xy$cell),colnames(obj)),max(abs(xy$x-input$x))==0,max(abs(xy$y-input$y))==0)
  obj <- NormalizeSMData(obj,normalisation.type="TIC",assay="SPM",verbose=FALSE)
  norm <- as.matrix(LayerData(obj,assay="SPM",layer="data"))
  expectation <- sweep(counts,2,colSums(counts),"/")*nrow(counts)
  err <- max(abs(norm-expectation))
  stopifnot(err<1e-8,all(is.finite(norm)))
  graph_edges <- NA_integer_;max_edge<-NA_real_;components<-NA_integer_
  if (binned) {
    featurevar <- apply(norm,1,var);features<-rownames(norm)[featurevar>0]
    obj <- ScaleData(obj,assay="SPM",features=features,verbose=FALSE)
    obj <- RunSpatialGraphPCA(obj,n_components=20,assay="SPM",slot="scale.data",image="fov",platform="ST",lambda=0.5,n_neighbors=NULL,verbose=FALSE)
    obj <- GetKmeanClusters(obj,reduction="SpatialPCA",cluster.name="SM_technical_region4",clusters=4L,iter.max=100L,nstart=20L,algorithm="Hartigan-Wong",seed=262235)
    emb<-Embeddings(obj,"SpatialPCA")
    stopifnot(nrow(emb)==ncol(obj),ncol(emb)==20,all(is.finite(emb)))
    g<-as(obj@graphs$SpatialKNN,"dgCMatrix")
    stopifnot(max(abs(g-t(g)))==0,all(diag(g)==0))
    graph_edges<-length(g@x)
    ij<-summary(g);max_edge<-max(sqrt((xy$x[ij$i]-xy$x[ij$j])^2+(xy$y[ij$i]-xy$y[ij$j])^2))
    components<-igraph::components(igraph::graph_from_adjacency_matrix(g,mode="undirected",weighted=TRUE))$no
    dd<-data.frame(cell=colnames(obj),x=xy$x,y=xy$y,region=obj$SM_technical_region4,emb,check.names=FALSE)
    write.table(dd,file.path(P,"tables",paste0(label,"_SPAMTP_GRAPH_REGIONS.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
    summary<-aggregate(t(norm),by=list(region=obj$SM_technical_region4),FUN=mean)
    write.table(summary,file.path(P,"tables",paste0(label,"_SPAMTP_REGION_TIC_MEANS.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
  }
  out<-file.path(P,"results","SpaMTP_objects",paste0(label,if(binned)"_bin4"else"_pixels","_SpaMTP.rds"))
  stopifnot(!file.exists(out))
  saveRDS(obj,out,compress="gzip")
  data.frame(label=label,binned=binned,pixels_or_bins=ncol(obj),features=nrow(obj),normalization_max_error=err,graph_edges=graph_edges,max_neighbor_distance_grid=max_edge,graph_components=components,rds_path=out,official_function_execution=TRUE,full_package_installed=FALSE)
}
for(label in c("R114-T-neg","R114-T-pos","Y12-T-neg","Y12-T-pos")){
 allqc[[length(allqc)+1]]<-run_one(label,FALSE);gc()
 allqc[[length(allqc)+1]]<-run_one(label,TRUE);gc()
}
qc<-do.call(rbind,allqc);print(qc,row.names=FALSE)
write.table(qc,file.path(P,"tables","SPAMTP_MODULE_QC.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(P,"environment","SPAMTP_MODULE_SESSION_INFO.txt"))
write_json(list(status="PASS",mode="unmodified_official_source_modules",objects=nrow(qc),no_new_figures=TRUE),file.path(P,"SPAMTP_MODULE_VALIDATION.json"),pretty=TRUE,auto_unbox=TRUE)
