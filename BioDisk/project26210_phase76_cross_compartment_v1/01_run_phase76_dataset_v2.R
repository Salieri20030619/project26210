#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 4L) stop('Usage: 01_run_phase76_dataset_v1.R <dataset> <seurat_rds> <metadata_tsv_gz> <outdir>')
dataset <- args[[1L]]
rds_path <- normalizePath(args[[2L]], winslash='/', mustWork=TRUE)
metadata_path <- normalizePath(args[[3L]], winslash='/', mustWork=TRUE)
outdir <- normalizePath(args[[4L]], winslash='/', mustWork=FALSE)
part <- paste0(outdir, '.part')
if (!(dataset %chin% c('GSE207493','GSE264586','GSE294109'))) stop('Unexpected dataset')
if (dir.exists(outdir) || file.exists(outdir) || dir.exists(part) || file.exists(part)) stop('Refusing existing output')
dir.create(part, recursive=TRUE)

TUMOR_MIN <- 30L
MYELOID_MIN <- 50L
FORMAL_MIN_PATIENTS <- 8L

PROGRAMS <- list(
  TUMOR_MTA_PRODUCTION_CONTEXT=c('MAT2A','MAT2B','AMD1','SRM','SMS'),
  MYELOID_MTA_RESPONSE=c('ARG1','VEGFA','IL10','TIMP1'),
  MYELOID_SPP1_STATE_NO_SPP1=c('GPNMB','LGALS3','FABP5','CD9','LPL'),
  MYELOID_IMMUNOSUPPRESSION_NO_IL10=c('TGFB1','VSIG4','MARCO','CD274','IDO1'),
  MYELOID_C1QC_STATE=c('C1QA','C1QB','C1QC'),
  MYELOID_MONOCYTE_PROGRAM=c('FCN1','VCAN','S100A8','S100A9','CTSS','LILRB1')
)
MIN_PRESENT <- c(TUMOR_MTA_PRODUCTION_CONTEXT=3L, MYELOID_MTA_RESPONSE=2L,
                 MYELOID_SPP1_STATE_NO_SPP1=3L, MYELOID_IMMUNOSUPPRESSION_NO_IL10=3L,
                 MYELOID_C1QC_STATE=2L, MYELOID_MONOCYTE_PROGRAM=4L)

RELATIONS <- data.table(
  relation_id=sprintf('P%02d',1:6),
  predictor=c('tumor_tymp_z_mean','myeloid_mta_response_mean','tumor_tymp_z_mean',
              'tumor_mta_production_mean','tumor_tymp_z_mean','tumor_tymp_z_mean'),
  outcome=c('myeloid_mta_response_mean','myeloid_spp1_state_mean','myeloid_spp1_state_mean',
            'myeloid_mta_response_mean','myeloid_spp1_detect_fraction','myeloid_immunosuppression_mean'),
  expected_direction=1L
)

atomic_fwrite <- function(x, path) {
  temp <- paste0(path,'.part')
  fwrite(x,temp,sep='\t',quote=TRUE,na='NA')
  if (!file.rename(temp,path)) stop('Atomic rename failed: ',path)
}

score_program <- function(expression, genes, min_present) {
  present <- intersect(genes,rownames(expression))
  if (length(present)<min_present) stop('Insufficient program coverage: ',paste(genes,collapse=';'))
  mat <- as.matrix(expression[present,,drop=FALSE])
  means <- rowMeans(mat)
  sds <- apply(mat,1,sd)
  variable <- is.finite(sds) & sds>0
  z <- matrix(0,nrow=nrow(mat),ncol=ncol(mat),dimnames=dimnames(mat))
  z[variable,] <- sweep(sweep(mat[variable,,drop=FALSE],1,means[variable],'-'),1,sds[variable],'/')
  list(score=colMeans(z),present=present,variable=present[variable])
}

gene_values <- function(expression,gene) {
  if (!(gene %chin% rownames(expression))) stop('Missing gene: ',gene)
  as.numeric(expression[gene,])
}

safe_cor <- function(x,y) {
  ok <- is.finite(x)&is.finite(y)
  if (sum(ok)<3L || uniqueN(x[ok])<3L || uniqueN(y[ok])<3L) return(c(rho=NA_real_,p=NA_real_))
  test <- suppressWarnings(cor.test(x[ok],y[ok],method='spearman',exact=FALSE))
  c(rho=unname(test$estimate),p=test$p.value)
}

loo_range <- function(x,y) {
  if (length(x)<4L) return(c(min=NA_real_,max=NA_real_))
  vals <- sapply(seq_along(x),function(i) suppressWarnings(cor(x[-i],y[-i],method='spearman',use='complete.obs')))
  c(min=min(vals,na.rm=TRUE),max=max(vals,na.rm=TRUE))
}

residual_cor <- function(x,y,n_tumor,n_myeloid) {
  ok <- is.finite(x)&is.finite(y)&is.finite(n_tumor)&is.finite(n_myeloid)
  if (sum(ok)<8L) return(c(rho=NA_real_,p=NA_real_))
  d <- data.frame(x=x[ok],y=y[ok],nt=log1p(n_tumor[ok]),nm=log1p(n_myeloid[ok]))
  xr <- residuals(lm(rank(x)~rank(nt)+rank(nm),data=d))
  yr <- residuals(lm(rank(y)~rank(nt)+rank(nm),data=d))
  test <- suppressWarnings(cor.test(xr,yr,method='spearman',exact=FALSE))
  c(rho=unname(test$estimate),p=test$p.value)
}

meta <- fread(metadata_path)
if (dataset=='GSE207493') {
  if (anyDuplicated(meta$cell_barcode)) stop('Duplicate metadata ID')
  meta[,`:=`(
    cell_id=cell_barcode,
    seurat_cell_name=paste0('GSE207493__GSE207493__',sub('^([^_]+)_','\\1__',cell_barcode)),
    patient_key=paste0('GSE207493::',sample),
    timepoint='Single',treatment='Unknown',response='Unknown',
    is_tumor=cell_type=='Tumor',
    is_broad_epithelial=cell_type=='Tumor',
    is_myeloid=cell_type %chin% c('Myeloid','Dendritic')
  )]
} else {
  if (anyDuplicated(meta$cell_id) || anyDuplicated(meta$seurat_cell_name)) stop('Duplicate metadata ID')
  meta[,`:=`(
    is_tumor=phase7_fine_annotation=='Epithelial_tumor_like',
    is_broad_epithelial=phase7_broad_annotation=='Epithelial',
    is_myeloid=phase7_broad_annotation=='Myeloid'
  )]
}
meta[,unit_id:=paste(patient_key,timepoint,sep='::')]

object <- readRDS(rds_path)
if (!inherits(object,'Seurat') || !('RNA' %chin% Assays(object))) stop('Invalid Seurat object')
eall <- LayerData(object[['RNA']],layer='data')
missing <- setdiff(meta$seurat_cell_name,colnames(eall))
if (length(missing)) stop('Metadata/Seurat mismatch: ',length(missing))

source_names <- rownames(eall); upper <- toupper(source_names)
collisions <- unique(upper[duplicated(upper)|duplicated(upper,fromLast=TRUE)])
needed <- unique(c('TYMP','MTAP','SPP1',unlist(PROGRAMS)))
available_idx <- which(upper %chin% setdiff(needed,collisions))
expression <- eall[available_idx,meta$seurat_cell_name,drop=FALSE]
rownames(expression) <- upper[available_idx]
colnames(expression) <- meta$cell_id
rm(object,eall); invisible(gc())

tumor_meta <- meta[is_tumor==TRUE]
broad_meta <- meta[is_broad_epithelial==TRUE]
myeloid_meta <- meta[is_myeloid==TRUE]
tumor_expr <- expression[,tumor_meta$cell_id,drop=FALSE]
broad_expr <- expression[,broad_meta$cell_id,drop=FALSE]
myeloid_expr <- expression[,myeloid_meta$cell_id,drop=FALSE]

tumor_tymp <- gene_values(tumor_expr,'TYMP')
tumor_tymp_sd <- sd(tumor_tymp)
tumor_tymp_z <- if (is.finite(tumor_tymp_sd) && tumor_tymp_sd>0) (tumor_tymp-mean(tumor_tymp))/tumor_tymp_sd else rep(0,length(tumor_tymp))
tumor_prod <- score_program(tumor_expr,PROGRAMS$TUMOR_MTA_PRODUCTION_CONTEXT,MIN_PRESENT[['TUMOR_MTA_PRODUCTION_CONTEXT']])

broad_tymp <- gene_values(broad_expr,'TYMP')
broad_tymp_sd <- sd(broad_tymp)
broad_tymp_z <- if (is.finite(broad_tymp_sd) && broad_tymp_sd>0) (broad_tymp-mean(broad_tymp))/broad_tymp_sd else rep(0,length(broad_tymp))

myeloid_scores <- list()
for (name in names(PROGRAMS)[names(PROGRAMS)!='TUMOR_MTA_PRODUCTION_CONTEXT']) {
  myeloid_scores[[name]] <- score_program(myeloid_expr,PROGRAMS[[name]],MIN_PRESENT[[name]])
}
myeloid_spp1 <- gene_values(myeloid_expr,'SPP1')

tumor_cells <- data.table(cell_id=tumor_meta$cell_id,unit_id=tumor_meta$unit_id,
                          tymp=tumor_tymp,tymp_z=tumor_tymp_z,mta_prod=tumor_prod$score)
broad_cells <- data.table(cell_id=broad_meta$cell_id,unit_id=broad_meta$unit_id,
                          tymp=broad_tymp,tymp_z=broad_tymp_z)
myeloid_cells <- data.table(cell_id=myeloid_meta$cell_id,unit_id=myeloid_meta$unit_id,
  mta_response=myeloid_scores$MYELOID_MTA_RESPONSE$score,
  spp1_state=myeloid_scores$MYELOID_SPP1_STATE_NO_SPP1$score,
  immunosuppression=myeloid_scores$MYELOID_IMMUNOSUPPRESSION_NO_IL10$score,
  c1qc_state=myeloid_scores$MYELOID_C1QC_STATE$score,
  monocyte_program=myeloid_scores$MYELOID_MONOCYTE_PROGRAM$score,
  spp1=myeloid_spp1)

base_units <- unique(meta[,.(unit_id,patient_key,timepoint,treatment,response)])
tc <- tumor_cells[,.(n_tumor_cells=.N,tumor_tymp_mean=mean(tymp),tumor_tymp_detect_fraction=mean(tymp>0),
                      tumor_tymp_z_mean=mean(tymp_z),tumor_mta_production_mean=mean(mta_prod)),by=unit_id]
bc <- broad_cells[,.(n_broad_epithelial_cells=.N,broad_epithelial_tymp_mean=mean(tymp),
                      broad_epithelial_tymp_detect_fraction=mean(tymp>0),broad_epithelial_tymp_z_mean=mean(tymp_z)),by=unit_id]
mc <- myeloid_cells[,.(n_myeloid_cells=.N,myeloid_mta_response_mean=mean(mta_response),
                        myeloid_spp1_state_mean=mean(spp1_state),myeloid_immunosuppression_mean=mean(immunosuppression),
                        myeloid_c1qc_state_mean=mean(c1qc_state),myeloid_monocyte_program_mean=mean(monocyte_program),
                        myeloid_spp1_mean=mean(spp1),myeloid_spp1_detect_fraction=mean(spp1>0)),by=unit_id]
unit <- Reduce(function(x,y) merge(x,y,by='unit_id',all.x=TRUE),list(base_units,tc,bc,mc))
for (col in c('n_tumor_cells','n_broad_epithelial_cells','n_myeloid_cells')) set(unit,which(is.na(unit[[col]])),col,0L)
unit[,eligible_primary:=n_tumor_cells>=TUMOR_MIN & n_myeloid_cells>=MYELOID_MIN]

score_cols <- c('tumor_tymp_mean','tumor_tymp_detect_fraction','tumor_tymp_z_mean','tumor_mta_production_mean',
                'broad_epithelial_tymp_mean','broad_epithelial_tymp_detect_fraction','broad_epithelial_tymp_z_mean',
                'myeloid_mta_response_mean','myeloid_spp1_state_mean','myeloid_immunosuppression_mean',
                'myeloid_c1qc_state_mean','myeloid_monocyte_program_mean','myeloid_spp1_mean','myeloid_spp1_detect_fraction')
patient <- unit[eligible_primary==TRUE, c(list(
  dataset=dataset,n_eligible_timepoints=.N,timepoints=paste(sort(unique(timepoint)),collapse=';'),
  treatments=paste(sort(unique(treatment)),collapse=';'),responses=paste(sort(unique(response)),collapse=';'),
  n_tumor_cells=sum(n_tumor_cells),n_myeloid_cells=sum(n_myeloid_cells)),
  lapply(.SD,function(x) median(x,na.rm=TRUE))),by=patient_key,.SDcols=score_cols]

primary_rows <- list(); sensitivity_rows <- list(); ri <- 0L; si <- 0L
for (i in seq_len(nrow(RELATIONS))) {
  rel <- RELATIONS[i]; x <- patient[[rel$predictor]]; y <- patient[[rel$outcome]]
  result <- safe_cor(x,y); loo <- loo_range(x,y); adjusted <- residual_cor(x,y,patient$n_tumor_cells,patient$n_myeloid_cells)
  ri <- ri+1L
  primary_rows[[ri]] <- data.table(dataset=dataset,analysis_family='PRIMARY_PATIENT',relation_id=rel$relation_id,
    predictor=rel$predictor,outcome=rel$outcome,expected_direction=1L,n_patients=nrow(patient),
    spearman_rho=result[['rho']],p_value=result[['p']],loo_min_rho=loo[['min']],loo_max_rho=loo[['max']],
    count_adjusted_rho=adjusted[['rho']],count_adjusted_p=adjusted[['p']])
}
primary <- rbindlist(primary_rows)
primary[,q_value:=p.adjust(p_value,method='BH')]
primary[,status:=fifelse(n_patients<FORMAL_MIN_PATIENTS,'DESCRIPTIVE_LOW_N',
  fifelse(q_value<0.05 & spearman_rho>0 & loo_min_rho>0,'COHORT_SUPPORTED',
    fifelse(q_value<0.05 & spearman_rho<=0,'DIRECTION_CONFLICT','NOT_SUPPORTED')))]

SENS <- data.table(
  sensitivity_id=c('S01','S02','S03','S04','S05','S06'),
  predictor=c('tumor_tymp_mean','tumor_tymp_detect_fraction','broad_epithelial_tymp_z_mean',
              'tumor_tymp_mean','tumor_tymp_detect_fraction','broad_epithelial_tymp_z_mean'),
  outcome=c('myeloid_mta_response_mean','myeloid_mta_response_mean','myeloid_mta_response_mean',
            'myeloid_spp1_state_mean','myeloid_spp1_state_mean','myeloid_spp1_state_mean'))
for (i in seq_len(nrow(SENS))) {
  rel <- SENS[i]; result <- safe_cor(patient[[rel$predictor]],patient[[rel$outcome]])
  si <- si+1L
  sensitivity_rows[[si]] <- data.table(dataset=dataset,analysis_family='SENSITIVITY_PATIENT',
    relation_id=rel$sensitivity_id,predictor=rel$predictor,outcome=rel$outcome,n_patients=nrow(patient),
    spearman_rho=result[['rho']],p_value=result[['p']])
}
sensitivity <- rbindlist(sensitivity_rows); sensitivity[,q_value:=p.adjust(p_value,method='BH')]

unit_rel <- list(); ui <- 0L
for (i in seq_len(nrow(RELATIONS))) {
  rel <- RELATIONS[i]; d <- unit[eligible_primary==TRUE]; result <- safe_cor(d[[rel$predictor]],d[[rel$outcome]])
  ui <- ui+1L; unit_rel[[ui]] <- data.table(dataset=dataset,analysis_family='DESCRIPTIVE_PATIENT_TIMEPOINT',
    relation_id=rel$relation_id,n_units=nrow(d),spearman_rho=result[['rho']],p_value=result[['p']])
}
unit_rel <- rbindlist(unit_rel); unit_rel[,q_value:=p.adjust(p_value,method='BH')]

coverage <- rbindlist(list(
  data.table(dataset=dataset,compartment='Tumor',program='TUMOR_MTA_PRODUCTION_CONTEXT',
             present_genes=paste(tumor_prod$present,collapse=';'),variable_genes=paste(tumor_prod$variable,collapse=';')),
  rbindlist(lapply(names(myeloid_scores),function(name) data.table(dataset=dataset,compartment='Myeloid',program=name,
    present_genes=paste(myeloid_scores[[name]]$present,collapse=';'),variable_genes=paste(myeloid_scores[[name]]$variable,collapse=';'))))
))

status <- data.table(dataset=dataset,source_rds=rds_path,source_metadata=metadata_path,
  source_object_cells=ncol(expression),metadata_cells=nrow(meta),tumor_cells=nrow(tumor_meta),
  broad_epithelial_cells=nrow(broad_meta),myeloid_cells=nrow(myeloid_meta),eligible_units=sum(unit$eligible_primary),
  eligible_patients=nrow(patient),formal_inference=nrow(patient)>=FORMAL_MIN_PATIENTS,
  uppercase_collision_keys=paste(collisions,collapse=';'),status='PASS_PHASE76_DATASET')

atomic_fwrite(unit,file.path(part,paste0(dataset,'_unit_summary_v1.tsv')))
atomic_fwrite(patient,file.path(part,paste0(dataset,'_patient_summary_v1.tsv')))
atomic_fwrite(primary,file.path(part,paste0(dataset,'_primary_relations_v1.tsv')))
atomic_fwrite(sensitivity,file.path(part,paste0(dataset,'_sensitivity_relations_v1.tsv')))
atomic_fwrite(unit_rel,file.path(part,paste0(dataset,'_unit_relations_descriptive_v1.tsv')))
atomic_fwrite(coverage,file.path(part,paste0(dataset,'_program_coverage_v1.tsv')))
atomic_fwrite(status,file.path(part,paste0(dataset,'_status_v1.tsv')))
if (!file.rename(part,outdir)) stop('Final atomic rename failed')
cat(sprintf('{"status":"PASS_PHASE76_DATASET","dataset":"%s","eligible_patients":%d}\n',dataset,nrow(patient)))

