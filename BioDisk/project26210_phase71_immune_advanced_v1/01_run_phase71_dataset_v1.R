#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) stop('Usage: 01_run_phase71_dataset_v1.R <dataset> <seurat_rds> <metadata_tsv_gz> <cycle_tsv> <output_dir>')
dataset <- args[[1L]]
seurat_path <- normalizePath(args[[2L]], winslash='/', mustWork=TRUE)
metadata_path <- normalizePath(args[[3L]], winslash='/', mustWork=TRUE)
cycle_path <- normalizePath(args[[4L]], winslash='/', mustWork=TRUE)
output_dir <- normalizePath(args[[5L]], winslash='/', mustWork=FALSE)
part_dir <- paste0(output_dir, '.part')
if (!dataset %chin% c('GSE207493','GSE264586','GSE294109')) stop('Unexpected dataset')
if (dir.exists(output_dir) || file.exists(output_dir) || dir.exists(part_dir) || file.exists(part_dir)) stop('Refusing existing output')
dir.create(part_dir, recursive=TRUE)

SEED <- 2621071L
MIN_CELLS <- 50L

PROGRAMS <- list(
  MTA_A2B_STAT3_RESPONSE = c('ARG1','VEGFA','IL10','TIMP1'),
  SPP1_STATE_NO_SPP1 = c('GPNMB','LGALS3','FABP5','CD9','LPL'),
  SPP1_STATE_FULL = c('SPP1','GPNMB','LGALS3','FABP5','CD9','LPL'),
  C1QC_STATE = c('C1QA','C1QB','C1QC'),
  MONOCYTE_PROGRAM = c('FCN1','VCAN','S100A8','S100A9','CTSS','LILRB1'),
  IMMUNOSUPPRESSION_NO_IL10 = c('TGFB1','VSIG4','MARCO','CD274','IDO1'),
  IMMUNOSUPPRESSION_WITH_IL10 = c('IL10','TGFB1','VSIG4','MARCO','CD274','IDO1'),
  MTA_PRODUCTION_CONTEXT = c('MAT2A','MAT2B','AMD1','SRM','SMS')
)
MIN_PRESENT <- c(MTA_A2B_STAT3_RESPONSE=2L, SPP1_STATE_NO_SPP1=3L, SPP1_STATE_FULL=4L,
                 C1QC_STATE=2L, MONOCYTE_PROGRAM=4L, IMMUNOSUPPRESSION_NO_IL10=3L,
                 IMMUNOSUPPRESSION_WITH_IL10=4L, MTA_PRODUCTION_CONTEXT=3L)

RELATIONS <- data.table(
  relation_id=c('P01','P02','P03','P04','P05','P06','S01','S02','S03','S04','S05','S06','S07'),
  tier=c(rep('PRIMARY',6),rep('SECONDARY',7)),
  predictor=c('TYMP','MTA_A2B_STAT3_RESPONSE','TYMP','MTA_PRODUCTION_CONTEXT','MTA_A2B_STAT3_RESPONSE','MTA_A2B_STAT3_RESPONSE','ADORA2B','STAT3','MTA_A2B_STAT3_RESPONSE','MTAP','MTA_A2B_STAT3_RESPONSE','MONOCYTE_PROGRAM','MTA_A2B_STAT3_RESPONSE'),
  outcome=c('MTA_A2B_STAT3_RESPONSE','SPP1_STATE_NO_SPP1','SPP1_STATE_NO_SPP1','MTA_A2B_STAT3_RESPONSE','IMMUNOSUPPRESSION_NO_IL10','C1QC_STATE','MTA_A2B_STAT3_RESPONSE','MTA_A2B_STAT3_RESPONSE','SPP1','MTA_A2B_STAT3_RESPONSE','SPP1_VS_C1QC','MTA_A2B_STAT3_RESPONSE','IMMUNOSUPPRESSION_WITH_IL10'),
  expected_direction=c('positive','positive','positive','positive','positive','two_sided','positive','positive','positive','negative','positive','negative','positive')
)

atomic_fwrite <- function(x, path) {
  temporary <- if (grepl('\\.gz$',path)) sub('\\.gz$','.part.gz',path) else paste0(path,'.part')
  fwrite(x, temporary, sep='\t', quote=TRUE, na='NA')
  if (!file.rename(temporary,path)) stop('Atomic rename failed: ',path)
}

score_matched_control <- function(expression, genes, ctrl_size, n_bins=25L, seed=SEED) {
  genes <- unique(intersect(toupper(genes), rownames(expression)))
  if (length(genes) < 10L) return(NULL)
  averages <- Matrix::rowMeans(expression)
  finite <- is.finite(averages)
  pool_names <- rownames(expression)[finite]
  pool_average <- averages[finite]
  n_items <- as.integer(round(length(pool_average)/(n_bins-1L)))
  gene_rank <- data.table::frank(pool_average,ties.method='min')
  bins <- floor(gene_rank/n_items); names(bins) <- pool_names
  set.seed(seed); controls <- character()
  for (bin in unique(bins[genes])) {
    candidates <- names(bins)[bins==bin & !names(bins) %chin% genes]
    if (length(candidates)>ctrl_size) candidates <- sample(candidates,ctrl_size,FALSE)
    controls <- union(controls,candidates)
  }
  if (!length(controls)) stop('No matched control genes')
  list(score=as.numeric(Matrix::colMeans(expression[genes,,drop=FALSE])-Matrix::colMeans(expression[controls,,drop=FALSE])), genes=genes, controls=controls)
}

residualize <- function(values,s_score,g2m_score) as.numeric(lm.fit(cbind(Intercept=1,S_score=s_score,G2M_score=g2m_score),values)$residuals)
safe_rho <- function(x,y) if (uniqueN(x)<3L || uniqueN(y)<3L) NA_real_ else suppressWarnings(cor(x,y,method='spearman'))
signed_rank <- function(values) {
  values <- values[is.finite(values)]
  if (!length(values) || all(abs(values)<1e-15)) return(c(statistic=0,p_value=1))
  z <- suppressWarnings(wilcox.test(values,mu=0,alternative='two.sided',exact=length(values)<50L))
  c(statistic=unname(z$statistic),p_value=z$p.value)
}

metadata <- fread(metadata_path)
object <- readRDS(seurat_path)
if (!inherits(object,'Seurat') || !'RNA' %chin% Assays(object)) stop('Invalid Seurat object')
expression_all <- LayerData(object[['RNA']],layer='data')
source_object_cells <- ncol(object)

if (dataset == 'GSE207493') {
  if (anyDuplicated(metadata$cell_barcode)) stop('Duplicate GSE207493 metadata IDs')
  metadata <- metadata[cell_type %chin% c('Myeloid','Dendritic')]
  metadata[, `:=`(cell_id=cell_barcode, seurat_cell_name=paste0('GSE207493__GSE207493__',sub('^([^_]+)_','\\1__',cell_barcode)),
                  patient_key=paste0('GSE207493::',sample), patient_id=sample,
                  timepoint='Single', treatment='Unknown', response='Unknown')]
} else {
  if (anyDuplicated(metadata$cell_id) || anyDuplicated(metadata$seurat_cell_name)) stop('Duplicate external metadata IDs')
  metadata <- metadata[phase7_broad_annotation=='Myeloid']
}
missing_cells <- setdiff(metadata$seurat_cell_name,colnames(expression_all))
if (length(missing_cells)) stop('Metadata/Seurat mismatch: ',length(missing_cells))
selected_columns <- match(metadata$seurat_cell_name,colnames(expression_all))
expression <- expression_all[,selected_columns,drop=FALSE]
colnames(expression) <- metadata$cell_id
rm(object,expression_all); invisible(gc())

source_names <- rownames(expression); upper_names <- toupper(source_names)
collision_keys <- unique(upper_names[duplicated(upper_names)|duplicated(upper_names,fromLast=TRUE)])
keep <- !upper_names %chin% collision_keys
expression <- expression[keep,,drop=FALSE]; rownames(expression) <- upper_names[keep]
if (anyDuplicated(rownames(expression)) || !'TYMP' %chin% rownames(expression)) stop('Invalid gene names')

cycle <- fread(cycle_path)
s_genes <- toupper(cycle[phase=='S',gene]); g2m_genes <- toupper(cycle[phase=='G2M',gene])
s_score <- score_matched_control(expression,s_genes,ctrl_size=min(length(s_genes),length(g2m_genes)))
g2m_score <- score_matched_control(expression,g2m_genes,ctrl_size=min(length(s_genes),length(g2m_genes)))
if (is.null(s_score)||is.null(g2m_score)) stop('Cycle scoring failed')

needed <- unique(c(unlist(PROGRAMS),RELATIONS$predictor,RELATIONS$outcome,'TYMP','MTAP','ADORA2B','STAT3','SPP1'))
present_needed <- intersect(needed,rownames(expression))
small <- as.matrix(expression[present_needed,,drop=FALSE])
gene_sd <- apply(small,1,sd)
z <- sweep(sweep(small,1,rowMeans(small),'-'),1,gene_sd,'/')
z[!is.finite(z)] <- 0
scores <- data.table(cell_id=colnames(expression))
coverage <- list(); idx <- 0L
for (name in names(PROGRAMS)) {
  genes <- intersect(PROGRAMS[[name]],rownames(z)); eligible <- length(genes)>=MIN_PRESENT[[name]]
  idx <- idx+1L; coverage[[idx]] <- data.table(dataset=dataset,program=name,n_frozen_genes=length(PROGRAMS[[name]]),n_present_genes=length(genes),present_genes=paste(genes,collapse=';'),eligible=eligible)
  if (!eligible) scores[[name]] <- NA_real_ else scores[[name]] <- colMeans(z[genes,,drop=FALSE])
}
for (gene in c('TYMP','MTAP','ADORA2B','STAT3','SPP1')) scores[[gene]] <- if (gene %chin% rownames(small)) as.numeric(small[gene,]) else NA_real_
scores[, S_score:=s_score$score]; scores[, G2M_score:=g2m_score$score]
scores[, SPP1_VS_C1QC:=SPP1_STATE_NO_SPP1-C1QC_STATE]
scores <- cbind(metadata[,.(cell_id,patient_key,patient_id,timepoint,treatment,response)],scores[,!"cell_id",with=FALSE])

rows <- list(); qrows <- list(); row_i <- 0L; q_i <- 0L
groups <- unique(scores[,.(patient_key,timepoint)])
for (g in seq_len(nrow(groups))) {
  key <- groups[g]; d <- scores[patient_key==key$patient_key & timepoint==key$timepoint]
  if (nrow(d)<MIN_CELLS) next
  for (r in seq_len(nrow(RELATIONS))) {
    rel <- RELATIONS[r]; x <- d[[rel$predictor]]; y <- d[[rel$outcome]]
    if (sum(is.finite(x)&is.finite(y))<MIN_CELLS || uniqueN(x[is.finite(x)])<3L || uniqueN(y[is.finite(y)])<3L) next
    ok <- is.finite(x)&is.finite(y)&is.finite(d$S_score)&is.finite(d$G2M_score)
    xr <- residualize(x[ok],d$S_score[ok],d$G2M_score[ok]); yr <- residualize(y[ok],d$S_score[ok],d$G2M_score[ok])
    row_i <- row_i+1L
    rows[[row_i]] <- data.table(dataset=dataset,relation_id=rel$relation_id,tier=rel$tier,predictor=rel$predictor,outcome=rel$outcome,expected_direction=rel$expected_direction,patient_key=key$patient_key,patient_id=d$patient_id[1],timepoint=key$timepoint,treatment=d$treatment[1],response=d$response[1],n_cells=sum(ok),raw_rho=safe_rho(x[ok],y[ok]),cycle_adjusted_rho=safe_rho(xr,yr))
    q25 <- quantile(x[ok],0.25,names=FALSE); q75 <- quantile(x[ok],0.75,names=FALSE)
    low <- y[ok][x[ok]<=q25]; high <- y[ok][x[ok]>=q75]
    if (length(low)>=10L && length(high)>=10L) {
      q_i <- q_i+1L; qrows[[q_i]] <- data.table(dataset=dataset,relation_id=rel$relation_id,tier=rel$tier,patient_key=key$patient_key,timepoint=key$timepoint,n_cells=sum(ok),n_low=length(low),n_high=length(high),top_minus_bottom=median(high)-median(low))
    }
  }
}
unit <- rbindlist(rows,use.names=TRUE,fill=TRUE)
quartile <- rbindlist(qrows,use.names=TRUE,fill=TRUE)
if (!nrow(unit)) stop('No eligible unit correlations')
patient <- unit[,.(patient_median_raw_rho=median(raw_rho,na.rm=TRUE),patient_median_adjusted_rho=median(cycle_adjusted_rho,na.rm=TRUE),n_timepoints=uniqueN(timepoint),total_cells=sum(n_cells)),by=.(dataset,relation_id,tier,predictor,outcome,expected_direction,patient_key)]
summary <- patient[,{
  test <- signed_rank(patient_median_adjusted_rho); pos <- sum(patient_median_adjusted_rho>0); neg <- sum(patient_median_adjusted_rho<0); nonzero <- pos+neg
  list(n_patients=uniqueN(patient_key),n_patient_timepoints=sum(n_timepoints),total_cells=sum(total_cells),median_raw_rho=median(patient_median_raw_rho),median_adjusted_rho=median(patient_median_adjusted_rho),positive_patients=pos,negative_patients=neg,positive_fraction=ifelse(nonzero>0,pos/nonzero,NA_real_),wilcoxon_statistic=test[['statistic']],p_value=test[['p_value']],sign_test_p=ifelse(nonzero>0,binom.test(pos,nonzero,p=0.5)$p.value,1))
},by=.(dataset,relation_id,tier,predictor,outcome,expected_direction)]
summary[,q_value:=p.adjust(p_value,method='BH'),by=.(dataset,tier)]
summary[,decision:=fifelse(expected_direction=='two_sided',fifelse(q_value<0.05,'ASSOCIATED_TWO_SIDED','NOT_SUPPORTED'),fifelse(q_value<0.05 & ((expected_direction=='positive' & median_adjusted_rho>0 & positive_fraction>=2/3)|(expected_direction=='negative' & median_adjusted_rho<0 & positive_fraction<=1/3)),'COHORT_SUPPORTED',fifelse(q_value<0.05,'DIRECTION_CONFLICT','NOT_SUPPORTED')))]

qpatient <- quartile[,.(patient_median_delta=median(top_minus_bottom,na.rm=TRUE),n_timepoints=uniqueN(timepoint)),by=.(dataset,relation_id,tier,patient_key)]
qsummary <- qpatient[,{test<-signed_rank(patient_median_delta); list(n_patients=uniqueN(patient_key),median_top_minus_bottom=median(patient_median_delta),p_value=test[['p_value']])},by=.(dataset,relation_id,tier)]
qsummary[,q_value:=p.adjust(p_value,method='BH'),by=.(dataset,tier)]

atomic_fwrite(rbindlist(coverage),file.path(part_dir,paste0(dataset,'_program_coverage_v1.tsv')))
atomic_fwrite(unit,file.path(part_dir,paste0(dataset,'_unit_correlations_v1.tsv.gz')))
atomic_fwrite(patient,file.path(part_dir,paste0(dataset,'_patient_correlations_v1.tsv')))
atomic_fwrite(summary,file.path(part_dir,paste0(dataset,'_relation_summary_v1.tsv')))
atomic_fwrite(quartile,file.path(part_dir,paste0(dataset,'_quartile_unit_effects_v1.tsv.gz')))
atomic_fwrite(qsummary,file.path(part_dir,paste0(dataset,'_quartile_summary_v1.tsv')))
atomic_fwrite(data.table(dataset=dataset,source_object=seurat_path,source_object_cells=source_object_cells,selected_cells=nrow(metadata),eligible_patients=uniqueN(patient$patient_key),genes_after_collision_exclusion=nrow(expression),uppercase_collision_keys=paste(collision_keys,collapse=';'),status='PASS_PHASE71_DATASET_V1'),file.path(part_dir,paste0(dataset,'_status_v1.tsv')))
file.rename(part_dir,output_dir)
cat(sprintf('{"status":"PASS_PHASE71_DATASET_V1","dataset":"%s","selected_cells":%d,"patients":%d}\n',dataset,nrow(metadata),uniqueN(patient$patient_key)))
