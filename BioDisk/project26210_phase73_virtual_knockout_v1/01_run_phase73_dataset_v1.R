#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
  library(scTenifoldKnk)
  library(scTenifoldNet)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4L) stop('Usage: 01_run_phase73_dataset_v1.R <dataset> <seurat_rds> <metadata_tsv_gz> <output_dir>')
dataset <- args[[1L]]
seurat_path <- normalizePath(args[[2L]],winslash='/',mustWork=TRUE)
metadata_path <- normalizePath(args[[3L]],winslash='/',mustWork=TRUE)
output_dir <- normalizePath(args[[4L]],winslash='/',mustWork=FALSE)
part_dir <- paste0(output_dir,'.part')

if (!dataset %in% c('GSE207493','GSE264586','GSE294109')) stop('Unexpected dataset')
if (dir.exists(output_dir) || file.exists(output_dir)) stop('Refusing existing completed output: ',output_dir)
dir.create(part_dir,recursive=TRUE,showWarnings=FALSE)

SEEDS <- c(GSE207493=26210731L,GSE264586=26210732L,GSE294109=26210733L)
SEED <- unname(SEEDS[[dataset]])
TARGETS <- c('STAT3','SPP1','TYMP','MTAP')
N_HVG <- 3000L
MIN_GENE_CELLS <- 25L
N_NET <- 10L
N_CELLS_NET <- 500L
N_COMP <- 3L
NETWORK_Q <- 0.90
TENSOR_K <- 3L
N_CORES <- 4L

PROGRAMS <- list(
  MTA_A2B_STAT3_RESPONSE=c('ARG1','VEGFA','IL10','TIMP1'),
  SPP1_STATE_NO_SPP1=c('GPNMB','LGALS3','FABP5','CD9','LPL'),
  IMMUNOSUPPRESSION_NO_IL10=c('TGFB1','VSIG4','MARCO','CD274','IDO1'),
  MTA_PRODUCTION_CONTEXT=c('MAT2A','MAT2B','AMD1','SRM','SMS'),
  C1QC_STATE=c('C1QA','C1QB','C1QC'),
  MONOCYTE_PROGRAM=c('FCN1','VCAN','S100A8','S100A9','CTSS','LILRB1'))
PRIMARY <- list(
  STAT3=c('MTA_A2B_STAT3_RESPONSE','SPP1_STATE_NO_SPP1','IMMUNOSUPPRESSION_NO_IL10'),
  TYMP=c('MTA_A2B_STAT3_RESPONSE','SPP1_STATE_NO_SPP1','IMMUNOSUPPRESSION_NO_IL10'),
  SPP1=c('SPP1_STATE_NO_SPP1','IMMUNOSUPPRESSION_NO_IL10'),
  MTAP=c('MTA_A2B_STAT3_RESPONSE','SPP1_STATE_NO_SPP1'))
FORCED <- unique(c(TARGETS,unlist(PROGRAMS,use.names=FALSE)))

atomic_fwrite <- function(x,path) {
  tmp <- if (grepl('\\.gz$',path)) sub('\\.gz$','.part.gz',path) else paste0(path,'.part')
  fwrite(x,tmp,sep='\t',quote=TRUE,na='NA')
  if (!file.rename(tmp,path)) stop('Atomic rename failed: ',path)
}
atomic_saveRDS <- function(x,path,compress=FALSE) {
  tmp <- paste0(path,'.part')
  saveRDS(x,tmp,compress=compress)
  if (!file.rename(tmp,path)) stop('Atomic rename failed: ',path)
}

message('LOAD_METADATA ',dataset)
metadata <- fread(metadata_path)
message('LOAD_SEURAT ',dataset)
object <- readRDS(seurat_path)
if (!inherits(object,'Seurat') || !'RNA' %in% Assays(object)) stop('Invalid Seurat object')
layers <- Layers(object[['RNA']])
count_layers <- layers[grepl('^counts($|\\.)',layers)]
if (length(count_layers) != 1L) stop('Expected exactly one counts layer')

if (dataset == 'GSE207493') {
  metadata <- metadata[cell_type %in% c('Myeloid','Dendritic')]
  metadata[,seurat_cell_name:=paste0('GSE207493__GSE207493__',sub('^([^_]+)_','\\1__',cell_barcode))]
} else {
  metadata <- metadata[phase7_broad_annotation=='Myeloid']
}
if (anyDuplicated(metadata$seurat_cell_name)) stop('Duplicate selected Seurat IDs')
missing <- setdiff(metadata$seurat_cell_name,colnames(object))
if (length(missing)) stop('Metadata/Seurat mismatch: ',length(missing))
counts <- LayerData(object[['RNA']],layer=count_layers[[1L]])[,metadata$seurat_cell_name,drop=FALSE]
source_object_cells <- ncol(object)
selected_cells_before_qc <- ncol(counts)
rm(object,metadata); invisible(gc())

upper <- toupper(rownames(counts))
collisions <- unique(upper[duplicated(upper)|duplicated(upper,fromLast=TRUE)])
keep <- !upper %in% collisions
counts <- counts[keep,,drop=FALSE]
rownames(counts) <- upper[keep]
if (anyDuplicated(rownames(counts))) stop('Uppercase gene collision remained')

message('SCQC ',dataset)
scqc <- get('scQC',envir=getNamespace('scTenifoldKnk'))
counts <- scqc(counts,mtThreshold=0.10,minLSize=1000)
gene_detected <- Matrix::rowSums(counts!=0)
counts <- counts[gene_detected>=MIN_GENE_CELLS,,drop=FALSE]
if (ncol(counts)<N_CELLS_NET) stop('Too few cells after QC')
if (!all(TARGETS %in% rownames(counts))) stop('Target missing after QC: ',paste(setdiff(TARGETS,rownames(counts)),collapse=';'))

message('SELECT_FEATURES ',dataset)
lib <- Matrix::colSums(counts)
norm <- counts %*% Diagonal(x=1/lib)
norm@x <- log1p(norm@x*1e4)
sq <- norm
sq@x <- sq@x^2
variance <- Matrix::rowMeans(sq)-Matrix::rowMeans(norm)^2
variance[!is.finite(variance)] <- -Inf
hvg <- names(sort(variance,decreasing=TRUE))[seq_len(min(N_HVG,length(variance)))]
model_genes <- unique(c(hvg,intersect(FORCED,rownames(counts))))
counts_model <- counts[model_genes,,drop=FALSE]
rm(norm,sq,counts); invisible(gc())

feature_table <- data.table(
  dataset=dataset,
  gene=rownames(counts_model),
  variance=variance[rownames(counts_model)],
  is_top_hvg=rownames(counts_model) %in% hvg,
  is_forced=rownames(counts_model) %in% FORCED,
  detected_cells=Matrix::rowSums(counts_model!=0))
atomic_fwrite(feature_table,file.path(part_dir,paste0(dataset,'_modeled_genes_v1.tsv.gz')))

wt_path <- file.path(part_dir,paste0(dataset,'_wt_network_v1.rds'))
if (file.exists(wt_path)) {
  message('LOAD_WT_CHECKPOINT ',dataset)
  WT <- as.matrix(readRDS(wt_path))
} else {
  message('BUILD_WT ',dataset,' genes=',nrow(counts_model),' cells=',ncol(counts_model))
  set.seed(SEED)
  nets <- scTenifoldNet::makeNetworks(
    X=counts_model,nNet=N_NET,nCells=N_CELLS_NET,nComp=N_COMP,
    scaleScores=TRUE,symmetric=FALSE,q=NETWORK_Q,nCores=N_CORES)
  message('TENSOR_DECOMPOSITION ',dataset)
  td <- scTenifoldNet::tensorDecomposition(
    xList=nets,K=TENSOR_K,maxError=1e-5,maxIter=1000,nDecimal=3)
  WT <- td$X
  strict_direction <- get('strictDirection',envir=getNamespace('scTenifoldKnk'))
  WT <- strict_direction(WT,lambda=0)
  WT <- as.matrix(WT)
  diag(WT) <- 0
  WT <- t(WT)
  atomic_saveRDS(Matrix(WT,sparse=TRUE),wt_path,compress=FALSE)
  rm(nets,td); invisible(gc())
}

dreg <- get('dRegulation',envir=getNamespace('scTenifoldKnk'))
program_files <- character()
for (target in TARGETS) {
  target_dir <- file.path(part_dir,paste0('KO_',target))
  dir.create(target_dir,recursive=TRUE,showWarnings=FALSE)
  marker <- file.path(target_dir,'KO_COMPLETE')
  program_path <- file.path(target_dir,paste0(dataset,'_',target,'_program_impact_v1.tsv'))
  diff_path <- file.path(target_dir,paste0(dataset,'_',target,'_diff_regulation_v1.tsv.gz'))
  if (file.exists(marker) && file.exists(program_path) && file.exists(diff_path)) {
    message('SKIP_COMPLETED_KO ',dataset,' ',target)
    program_files <- c(program_files,program_path)
    next
  }
  message('SIMULATE_KO ',dataset,' ',target)
  KO <- WT
  KO[target,] <- 0
  MA <- scTenifoldNet::manifoldAlignment(WT,KO,d=2,nCores=N_CORES)
  DR <- as.data.table(dreg(MA,target))
  setnames(DR,'gene','affected_gene')
  DR[,`:=`(dataset=dataset,knockout_gene=target)]
  DR[,impact_rank:=frank(-distance,ties.method='average')]
  DR[,impact_percentile:=frank(distance,ties.method='average')/.N]
  setcolorder(DR,c('dataset','knockout_gene','affected_gene','distance','Z','FC','p.value','p.adj','impact_rank','impact_percentile'))
  atomic_fwrite(DR,diff_path)

  prows <- list()
  for (program in names(PROGRAMS)) {
    pg <- intersect(PROGRAMS[[program]],DR$affected_gene)
    pg <- setdiff(pg,target)
    bg <- DR[affected_gene!=target & !affected_gene %in% pg,distance]
    pv <- DR[affected_gene %in% pg,distance]
    pp <- DR[affected_gene %in% pg,impact_percentile]
    if (length(pg)>=2L && length(bg)>20L && length(unique(c(pv,bg)))>1L) {
      test <- suppressWarnings(wilcox.test(pv,bg,alternative='greater',exact=FALSE))
      pval <- test$p.value
      statistic <- unname(test$statistic)
    } else {
      pval <- NA_real_; statistic <- NA_real_
    }
    prows[[length(prows)+1L]] <- data.table(
      dataset=dataset,knockout_gene=target,program=program,
      hypothesis_tier=ifelse(program %in% PRIMARY[[target]],'PRIMARY','CONTEXT'),
      n_frozen_genes=length(PROGRAMS[[program]]),n_modeled_genes=length(pg),
      modeled_genes=paste(pg,collapse=';'),median_distance=ifelse(length(pv),median(pv),NA_real_),
      median_impact_percentile=ifelse(length(pp),median(pp),NA_real_),
      top10_fraction=ifelse(length(pp),mean(pp>=0.90),NA_real_),
      wilcoxon_statistic=statistic,p_value=pval)
  }
  P <- rbindlist(prows)
  P[,q_value:=p.adjust(p_value,method='BH')]
  P[,decision:=ifelse(is.finite(q_value) & q_value<0.05 & median_impact_percentile>=0.75,
                       'PROGRAM_IMPACT_SUPPORTED','NOT_SUPPORTED')]
  atomic_fwrite(P,program_path)
  writeLines(c('KO_COMPLETE',paste0('dataset=',dataset),paste0('target=',target),
               paste0('completed_at=',format(Sys.time(),'%Y-%m-%dT%H:%M:%S%z'))),marker)
  program_files <- c(program_files,program_path)
  rm(KO,MA,DR,P); invisible(gc())
}

program_all <- rbindlist(lapply(program_files,fread),use.names=TRUE,fill=TRUE)
atomic_fwrite(program_all,file.path(part_dir,paste0(dataset,'_program_impact_all_v1.tsv')))
status <- data.table(
  dataset=dataset,status='PASS_PHASE73_DATASET_V1',source_rds=seurat_path,
  metadata=metadata_path,source_object_cells=source_object_cells,
  selected_cells_before_qc=selected_cells_before_qc,retained_cells_after_qc=ncol(counts_model),
  modeled_genes=nrow(counts_model),targets=paste(TARGETS,collapse=';'),
  n_networks=N_NET,cells_per_network=N_CELLS_NET,network_q=NETWORK_Q,tensor_k=TENSOR_K,
  seed=SEED,r_version=R.version.string,scTenifoldKnk_version=as.character(packageVersion('scTenifoldKnk')))
atomic_fwrite(status,file.path(part_dir,paste0(dataset,'_status_v1.tsv')))
writeLines(c('DATASET_COMPLETE',paste0('dataset=',dataset)),file.path(part_dir,'DATASET_COMPLETE'))
if (!file.rename(part_dir,output_dir)) stop('Failed to finalize dataset directory')
cat(sprintf('PASS_PHASE73_DATASET_V1 dataset=%s cells=%d genes=%d\n',dataset,ncol(counts_model),nrow(counts_model)))
