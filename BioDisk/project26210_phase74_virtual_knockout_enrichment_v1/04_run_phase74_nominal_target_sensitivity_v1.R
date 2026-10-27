#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=3L) stop('Usage: 04_run_phase74_nominal_target_sensitivity_v1.R <phase73_dataset_root> <phase74_enrichment_all_tsv_gz> <output_dir>')
phase73_root <- normalizePath(args[[1L]],winslash='/',mustWork=TRUE)
enrichment_path <- normalizePath(args[[2L]],winslash='/',mustWork=TRUE)
output_dir <- normalizePath(args[[3L]],winslash='/',mustWork=FALSE)
part_dir <- paste0(output_dir,'.part')
if (dir.exists(output_dir)||file.exists(output_dir)||dir.exists(part_dir)||file.exists(part_dir)) stop('Refusing existing output')
dir.create(part_dir,recursive=TRUE)

TARGETS <- c('TYMP','MTAP','SPP1','STAT3')
DATASETS <- c('GSE207493','GSE264586','GSE294109')
KOS <- c('STAT3','SPP1','TYMP','MTAP')

atomic_fwrite <- function(x,path) {
  tmp <- paste0(path,'.part'); fwrite(x,tmp,sep='\t',quote=TRUE,na='NA')
  if (!file.rename(tmp,path)) stop('Atomic rename failed: ',path)
}

files <- list.files(phase73_root,pattern='_diff_regulation_v1\\.tsv\\.gz$',recursive=TRUE,full.names=TRUE)
if (length(files)!=12L) stop('Expected 12 differential-regulation files')
grows <- list()
for (f in sort(files)) {
  d <- fread(f); dataset <- unique(d$dataset); ko <- unique(d$knockout_gene)
  if (length(dataset)!=1L||length(ko)!=1L) stop('Mixed input')
  for (target in TARGETS) {
    x <- d[affected_gene==target]
    if (nrow(x)!=1L) stop('Expected one target row: ',dataset,' ',ko,' ',target)
    grows[[length(grows)+1L]] <- data.table(
      dataset=dataset,knockout_gene=ko,queried_target_gene=target,
      self_knockout=ko==target,distance=x$distance,Z=x$Z,FC=x$FC,
      p_value=x$p.value,q_value=x$p.adj,impact_rank=x$impact_rank,
      impact_percentile=x$impact_percentile,
      nominal_gene_significant=is.finite(x$p.value)&x$p.value<0.05,
      fdr_gene_significant=is.finite(x$p.adj)&x$p.adj<0.05)
  }
}
gene_audit <- rbindlist(grows)
atomic_fwrite(gene_audit,file.path(part_dir,'phase74_gene_level_target_audit_v1.tsv'))

gcross <- list()
for (ko in KOS) for (target in TARGETS) {
  x <- gene_audit[knockout_gene==ko&queried_target_gene==target]
  getrow <- function(ds) x[dataset==ds][1]
  d0<-getrow('GSE207493');d1<-getrow('GSE264586');d2<-getrow('GSE294109')
  ne <- sum(c(d1$nominal_gene_significant,d2$nominal_gene_significant))
  nf <- sum(c(d1$fdr_gene_significant,d2$fdr_gene_significant))
  gcross[[length(gcross)+1L]] <- data.table(
    knockout_gene=ko,queried_target_gene=target,self_knockout=ko==target,
    external_nominal_decision=ifelse(ne==2,'TWO_EXTERNAL_GENE_NOMINAL',ifelse(ne==1,'ONE_EXTERNAL_GENE_NOMINAL','NOT_EXTERNALLY_NOMINAL')),
    external_fdr_decision=ifelse(nf==2,'TWO_EXTERNAL_GENE_FDR',ifelse(nf==1,'ONE_EXTERNAL_GENE_FDR','NOT_EXTERNALLY_FDR')),
    GSE207493_p=d0$p_value,GSE207493_q=d0$q_value,GSE207493_nominal=d0$nominal_gene_significant,GSE207493_fdr=d0$fdr_gene_significant,
    GSE264586_p=d1$p_value,GSE264586_q=d1$q_value,GSE264586_nominal=d1$nominal_gene_significant,GSE264586_fdr=d1$fdr_gene_significant,
    GSE294109_p=d2$p_value,GSE294109_q=d2$q_value,GSE294109_nominal=d2$nominal_gene_significant,GSE294109_fdr=d2$fdr_gene_significant)
}
gene_cross <- rbindlist(gcross)
atomic_fwrite(gene_cross,file.path(part_dir,'phase74_gene_level_target_cross_cohort_v1.tsv'))

enrich <- fread(enrichment_path)
split_genes <- strsplit(enrich$affected_genes,'/',fixed=TRUE)
enrich[,matched_target_genes:=vapply(split_genes,function(x)paste(intersect(TARGETS,x),collapse=';'),character(1))]
enrich[,nominal_pathway_significant:=is.finite(p_value)&p_value<0.05&affected_gene_count>=3L]
enrich[,fdr_pathway_significant:=is.finite(q_value)&q_value<0.05&affected_gene_count>=3L]
target_rows <- enrich[nominal_pathway_significant==TRUE&nzchar(matched_target_genes)]
atomic_fwrite(target_rows,file.path(part_dir,'phase74_nominal_enrichment_target_gene_rows_v1.tsv'))

keys <- unique(target_rows[,.(knockout_gene,ontology,term_id,term_name)])
pcross <- list()
for (i in seq_len(nrow(keys))) {
  k <- keys[i]; x <- target_rows[knockout_gene==k$knockout_gene&ontology==k$ontology&term_id==k$term_id]
  getrow <- function(ds) x[dataset==ds][1]
  d0<-getrow('GSE207493');d1<-getrow('GSE264586');d2<-getrow('GSE294109')
  present <- c(nrow(d1)>0,nrow(d2)>0)
  g1 <- if(nrow(d1))strsplit(d1$matched_target_genes,';',fixed=TRUE)[[1]] else character()
  g2 <- if(nrow(d2))strsplit(d2$matched_target_genes,';',fixed=TRUE)[[1]] else character()
  common <- intersect(g1,g2)
  decision <- if (all(present)&&length(common)) 'TWO_EXTERNAL_NOMINAL_SAME_TARGET' else if (all(present)) 'TWO_EXTERNAL_NOMINAL_DIFFERENT_TARGET' else if (sum(present)==1) 'ONE_EXTERNAL_NOMINAL_TARGET' else 'DISCOVERY_ONLY_NOMINAL_TARGET'
  pcross[[length(pcross)+1L]] <- data.table(
    knockout_gene=k$knockout_gene,ontology=k$ontology,term_id=k$term_id,term_name=k$term_name,
    external_nominal_decision=decision,common_external_target_genes=paste(common,collapse=';'),
    GSE207493_p=if(nrow(d0))d0$p_value else NA_real_,GSE207493_q=if(nrow(d0))d0$q_value else NA_real_,GSE207493_targets=if(nrow(d0))d0$matched_target_genes else '',
    GSE264586_p=if(nrow(d1))d1$p_value else NA_real_,GSE264586_q=if(nrow(d1))d1$q_value else NA_real_,GSE264586_targets=if(nrow(d1))d1$matched_target_genes else '',
    GSE294109_p=if(nrow(d2))d2$p_value else NA_real_,GSE294109_q=if(nrow(d2))d2$q_value else NA_real_,GSE294109_targets=if(nrow(d2))d2$matched_target_genes else '')
}
path_cross <- rbindlist(pcross,use.names=TRUE,fill=TRUE)
setorder(path_cross,knockout_gene,ontology,external_nominal_decision,term_id)
atomic_fwrite(path_cross,file.path(part_dir,'phase74_nominal_enrichment_target_gene_cross_cohort_v1.tsv'))

term_counts <- rbindlist(lapply(KOS,function(ko)rbindlist(lapply(TARGETS,function(target){
  rows <- target_rows[knockout_gene==ko&vapply(strsplit(matched_target_genes,';',fixed=TRUE),function(x)target%in%x,logical(1))]
  cross <- path_cross[knockout_gene==ko&vapply(strsplit(common_external_target_genes,';',fixed=TRUE),function(x)target%in%x,logical(1))]
  data.table(knockout_gene=ko,target_gene_in_overlap=target,
             GSE207493_nominal_terms=sum(rows$dataset=='GSE207493'),GSE264586_nominal_terms=sum(rows$dataset=='GSE264586'),GSE294109_nominal_terms=sum(rows$dataset=='GSE294109'),
             two_external_same_target_terms=sum(cross$external_nominal_decision=='TWO_EXTERNAL_NOMINAL_SAME_TARGET'))
}))))
atomic_fwrite(term_counts,file.path(part_dir,'phase74_nominal_target_gene_term_counts_v1.tsv'))

summary <- c('# Phase74 nominal target-gene sensitivity summary v1','',
             'Nominal P<0.05 is reported separately and does not replace BH-FDR decisions. All effects are unsigned network effects.','',
             '## Gene-level cross-cohort target audit','')
for (i in seq_len(nrow(gene_cross))) {
  x<-gene_cross[i]
  summary<-c(summary,sprintf('- %s KO -> %s: nominal=%s; FDR=%s; self_KO=%s.',x$knockout_gene,x$queried_target_gene,x$external_nominal_decision,x$external_fdr_decision,x$self_knockout))
}
summary<-c(summary,'','## Pathway overlap target-gene audit','',sprintf('- Nominal pathway rows containing at least one target gene: %d.',nrow(target_rows)),sprintf('- Cross-cohort unique target-containing terms: %d.',nrow(path_cross)),sprintf('- Two-external nominal terms containing the same target gene: %d.',sum(path_cross$external_nominal_decision=='TWO_EXTERNAL_NOMINAL_SAME_TARGET')),'')
writeLines(summary,file.path(part_dir,'PHASE74_NOMINAL_TARGET_GENE_SUMMARY_v1.md'))
metadata<-data.table(status='PASS_PHASE74_NOMINAL_TARGET_SENSITIVITY_V1',gene_audit_rows=nrow(gene_audit),gene_cross_rows=nrow(gene_cross),nominal_target_pathway_rows=nrow(target_rows),pathway_cross_rows=nrow(path_cross),two_external_same_target_terms=sum(path_cross$external_nominal_decision=='TWO_EXTERNAL_NOMINAL_SAME_TARGET'),direction_available=FALSE)
atomic_fwrite(metadata,file.path(part_dir,'run_metadata_v1.tsv'))
writeLines('PHASE74_NOMINAL_TARGET_SENSITIVITY_PASS',file.path(part_dir,'NOMINAL_SENSITIVITY_COMPLETE'))
if(!file.rename(part_dir,output_dir))stop('Failed to finalize nominal sensitivity')
cat(sprintf('PASS_PHASE74_NOMINAL_TARGET gene_rows=%d pathway_rows=%d cross=%d\n',nrow(gene_audit),nrow(target_rows),nrow(path_cross)))
