#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2L)stop('Usage: 05_finalize_phase74_nominal_target_sensitivity_v2.R <v1_dir> <v2_output_dir>')
v1_dir<-normalizePath(args[[1]],winslash='/',mustWork=TRUE)
out<-normalizePath(args[[2]],winslash='/',mustWork=FALSE);part<-paste0(out,'.part')
if(dir.exists(out)||file.exists(out)||dir.exists(part)||file.exists(part))stop('Refusing existing output')
dir.create(part,recursive=TRUE)

TARGETS<-c('TYMP','MTAP','SPP1','STAT3');KOS<-c('STAT3','SPP1','TYMP','MTAP');DATASETS<-c('GSE207493','GSE264586','GSE294109')
atomic_fwrite<-function(x,path){tmp<-paste0(path,'.part');fwrite(x,tmp,sep='\t',quote=TRUE,na='NA');if(!file.rename(tmp,path))stop('rename failed')}

gene_audit<-fread(file.path(v1_dir,'phase74_gene_level_target_audit_v1.tsv'))
gene_cross<-fread(file.path(v1_dir,'phase74_gene_level_target_cross_cohort_v1.tsv'))
target_rows<-fread(file.path(v1_dir,'phase74_nominal_enrichment_target_gene_rows_v1.tsv'))
if(nrow(gene_audit)!=48L||nrow(gene_cross)!=16L||nrow(target_rows)!=587L)stop('Unexpected v1 source dimensions')
atomic_fwrite(gene_audit,file.path(part,'phase74_gene_level_target_audit_v2.tsv'))
atomic_fwrite(gene_cross,file.path(part,'phase74_gene_level_target_cross_cohort_v2.tsv'))
atomic_fwrite(target_rows,file.path(part,'phase74_nominal_enrichment_target_gene_rows_v2.tsv'))

safe_row<-function(x,ds){y<-x[dataset==ds];if(!nrow(y))return(NULL);y[1]}
safe_targets<-function(x){if(is.null(x)||is.na(x$matched_target_genes)||!nzchar(x$matched_target_genes))character() else strsplit(x$matched_target_genes,';',fixed=TRUE)[[1]]}
keys<-unique(target_rows[,.(knockout_gene,ontology,term_id,term_name)])
rows<-list()
for(i in seq_len(nrow(keys))){
  k<-keys[i];x<-target_rows[knockout_gene==k$knockout_gene&ontology==k$ontology&term_id==k$term_id]
  d0<-safe_row(x,'GSE207493');d1<-safe_row(x,'GSE264586');d2<-safe_row(x,'GSE294109')
  p1<-!is.null(d1);p2<-!is.null(d2);g1<-safe_targets(d1);g2<-safe_targets(d2);common<-intersect(g1,g2)
  decision<-if(p1&&p2&&length(common))'TWO_EXTERNAL_NOMINAL_SAME_TARGET' else if(p1&&p2)'TWO_EXTERNAL_NOMINAL_DIFFERENT_TARGET' else if(xor(p1,p2))'ONE_EXTERNAL_NOMINAL_TARGET' else 'DISCOVERY_ONLY_NOMINAL_TARGET'
  rows[[length(rows)+1L]]<-data.table(
    knockout_gene=k$knockout_gene,ontology=k$ontology,term_id=k$term_id,term_name=k$term_name,
    external_nominal_decision=decision,common_external_target_genes=paste(common,collapse=';'),
    GSE207493_present=!is.null(d0),GSE207493_p=if(is.null(d0))NA_real_ else d0$p_value,GSE207493_q=if(is.null(d0))NA_real_ else d0$q_value,GSE207493_targets=paste(safe_targets(d0),collapse=';'),
    GSE264586_present=p1,GSE264586_p=if(is.null(d1))NA_real_ else d1$p_value,GSE264586_q=if(is.null(d1))NA_real_ else d1$q_value,GSE264586_targets=paste(g1,collapse=';'),
    GSE294109_present=p2,GSE294109_p=if(is.null(d2))NA_real_ else d2$p_value,GSE294109_q=if(is.null(d2))NA_real_ else d2$q_value,GSE294109_targets=paste(g2,collapse=';'))
}
cross<-rbindlist(rows);setorder(cross,knockout_gene,ontology,external_nominal_decision,term_id)
atomic_fwrite(cross,file.path(part,'phase74_nominal_enrichment_target_gene_cross_cohort_v2.tsv'))

contains_target<-function(strings,target)vapply(strsplit(strings,';',fixed=TRUE),function(x)target%in%x,logical(1))
counts<-rbindlist(lapply(KOS,function(ko)rbindlist(lapply(TARGETS,function(target){
  rr<-target_rows[knockout_gene==ko&contains_target(matched_target_genes,target)]
  cc<-cross[knockout_gene==ko&external_nominal_decision=='TWO_EXTERNAL_NOMINAL_SAME_TARGET'&contains_target(common_external_target_genes,target)]
  data.table(knockout_gene=ko,target_gene_in_overlap=target,GSE207493_nominal_terms=sum(rr$dataset=='GSE207493'),GSE264586_nominal_terms=sum(rr$dataset=='GSE264586'),GSE294109_nominal_terms=sum(rr$dataset=='GSE294109'),two_external_same_target_terms=nrow(cc))
}))))
atomic_fwrite(counts,file.path(part,'phase74_nominal_target_gene_term_counts_v2.tsv'))

nonself<-gene_cross[self_knockout==FALSE&(external_nominal_decision!='NOT_EXTERNALLY_NOMINAL'|GSE207493_nominal==TRUE)]
summary<-c('# Phase74 nominal target-gene sensitivity summary v2','',
  'Nominal P<0.05 is a sensitivity layer and does not replace BH-FDR. Effects are unsigned. Self-KO target hits are structural and not downstream evidence.','',
  '## Non-self gene-level nominal findings','')
if(nrow(nonself))for(i in seq_len(nrow(nonself))){x<-nonself[i];summary<-c(summary,sprintf('- %s KO -> %s: %s; external P values %.6g / %.6g; external q values %.6g / %.6g.',x$knockout_gene,x$queried_target_gene,x$external_nominal_decision,x$GSE264586_p,x$GSE294109_p,x$GSE264586_q,x$GSE294109_q))} else summary<-c(summary,'- None.')
summary<-c(summary,'','## Pathway target membership','',sprintf('- Nominal target-containing cohort rows: %d.',nrow(target_rows)),sprintf('- Unique cross-cohort terms: %d.',nrow(cross)),sprintf('- Two-external terms containing the same target: %d.',sum(cross$external_nominal_decision=='TWO_EXTERNAL_NOMINAL_SAME_TARGET')),'',
  'All target genes found in pathway overlap lists are reported in the count table; self-target membership must not be treated as independent mechanistic validation.','')
writeLines(summary,file.path(part,'PHASE74_NOMINAL_TARGET_GENE_SUMMARY_v2.md'))
metadata<-data.table(status='PASS_PHASE74_NOMINAL_TARGET_SENSITIVITY_V2_VALIDATION_PENDING',gene_audit_rows=nrow(gene_audit),gene_cross_rows=nrow(gene_cross),nominal_target_pathway_rows=nrow(target_rows),pathway_cross_rows=nrow(cross),two_external_same_target_terms=sum(cross$external_nominal_decision=='TWO_EXTERNAL_NOMINAL_SAME_TARGET'),direction_available=FALSE,v1_status='INVALID_CROSS_COHORT_MISSING_ROW_BUG')
atomic_fwrite(metadata,file.path(part,'run_metadata_v2.tsv'))
writeLines('PHASE74_NOMINAL_TARGET_SENSITIVITY_V2_PASS',file.path(part,'NOMINAL_SENSITIVITY_COMPLETE_V2'))
if(!file.rename(part,out))stop('Failed final rename')
cat(sprintf('PASS_PHASE74_NOMINAL_TARGET_V2 cross=%d two_external_same=%d\n',nrow(cross),sum(cross$external_nominal_decision=='TWO_EXTERNAL_NOMINAL_SAME_TARGET')))
