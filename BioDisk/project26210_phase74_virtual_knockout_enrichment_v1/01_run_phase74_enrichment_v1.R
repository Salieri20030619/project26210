#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(msigdbr)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop('Usage: 01_run_phase74_enrichment_v1.R <phase73_dataset_root> <output_dir>')
input_root <- normalizePath(args[[1L]],winslash='/',mustWork=TRUE)
output_dir <- normalizePath(args[[2L]],winslash='/',mustWork=FALSE)
part_dir <- paste0(output_dir,'.part')
if (dir.exists(output_dir)||file.exists(output_dir)||dir.exists(part_dir)||file.exists(part_dir)) stop('Refusing existing output')
dir.create(part_dir,recursive=TRUE)

atomic_fwrite <- function(x,path) {
  tmp <- if (grepl('\\.gz$',path)) sub('\\.gz$','.part.gz',path) else paste0(path,'.part')
  fwrite(x,tmp,sep='\t',quote=TRUE,na='NA')
  if (!file.rename(tmp,path)) stop('Atomic rename failed: ',path)
}

standardize_enrichment <- function(x,dataset,target,ontology) {
  if (is.null(x)) return(data.table())
  d <- as.data.table(as.data.frame(x))
  if (!nrow(d)) return(data.table())
  data.table(
    dataset=dataset,knockout_gene=target,ontology=ontology,
    term_id=as.character(d$ID),term_name=as.character(d$Description),
    gene_ratio=as.character(d$GeneRatio),background_ratio=as.character(d$BgRatio),
    p_value=as.numeric(d$pvalue),q_value=as.numeric(d$p.adjust),
    affected_gene_count=as.integer(d$Count),affected_genes=as.character(d$geneID))
}

message('LOAD_OFFLINE_KEGG_CACHE')
kegg <- as.data.table(msigdbr(db_species='HS',species='Homo sapiens',collection='C2'))
kegg <- kegg[gs_subcollection %in% c('CP:KEGG_LEGACY','CP:KEGG_MEDICUS')]
if (!nrow(kegg)) stop('No cached KEGG sets')
kegg[,ontology:=fifelse(gs_subcollection=='CP:KEGG_LEGACY','KEGG_LEGACY','KEGG_MEDICUS')]
kegg_sets <- split(kegg,kegg$ontology)

files <- list.files(input_root,pattern='_diff_regulation_v1\\.tsv\\.gz$',recursive=TRUE,full.names=TRUE)
if (length(files)!=12L) stop('Expected 12 Phase73 differential-regulation files, found ',length(files))

affected_all <- list(); count_rows <- list(); enrich_all <- list()
for (f in sort(files)) {
  d <- fread(f)
  if (!all(c('dataset','knockout_gene','affected_gene','distance','Z','FC','p.value','p.adj','impact_rank','impact_percentile') %in% names(d))) stop('Invalid input: ',f)
  dataset <- unique(d$dataset); target <- unique(d$knockout_gene)
  if (length(dataset)!=1L||length(target)!=1L) stop('Mixed dataset/target: ',f)
  universe_symbols <- unique(d$affected_gene)
  affected <- d[is.finite(p.adj)&p.adj<0.05]
  if (nrow(affected)<10L) stop('Too few affected genes: ',dataset,' ',target)
  affected[,source_file:=f]
  affected_all[[length(affected_all)+1L]] <- affected

  map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,keys=universe_symbols,columns='ENTREZID',keytype='SYMBOL'))
  map <- unique(as.data.table(map)[!is.na(ENTREZID),.(SYMBOL,ENTREZID=as.character(ENTREZID))])
  universe_entrez <- unique(map$ENTREZID)
  affected_entrez <- unique(map[SYMBOL %in% affected$affected_gene,ENTREZID])
  count_rows[[length(count_rows)+1L]] <- data.table(
    dataset=dataset,knockout_gene=target,modeled_genes=nrow(d),affected_fdr05=nrow(affected),
    affected_fraction=nrow(affected)/nrow(d),universe_entrez=length(universe_entrez),
    affected_entrez=length(affected_entrez),symbol_mapping_fraction=uniqueN(map$SYMBOL)/length(universe_symbols))

  message('ENRICH ',dataset,' ',target,' affected=',nrow(affected))
  for (ont in c('BP','CC','MF')) {
    e <- suppressMessages(enrichGO(
      gene=affected_entrez,universe=universe_entrez,OrgDb=org.Hs.eg.db,
      keyType='ENTREZID',ont=ont,pAdjustMethod='BH',pvalueCutoff=1,qvalueCutoff=1,
      minGSSize=10,maxGSSize=500,readable=TRUE))
    enrich_all[[length(enrich_all)+1L]] <- standardize_enrichment(e,dataset,target,paste0('GO_',ont))
  }
  for (ont in names(kegg_sets)) {
    ks <- kegg_sets[[ont]]
    term2gene <- unique(ks[,.(term=gs_name,gene=gene_symbol)])
    term2name <- unique(ks[,.(term=gs_name,name=gs_name)])
    e <- suppressMessages(enricher(
      gene=affected$affected_gene,universe=universe_symbols,
      TERM2GENE=term2gene,TERM2NAME=term2name,pAdjustMethod='BH',
      pvalueCutoff=1,qvalueCutoff=1,minGSSize=10,maxGSSize=500))
    enrich_all[[length(enrich_all)+1L]] <- standardize_enrichment(e,dataset,target,ont)
  }
}

affected_dt <- rbindlist(affected_all,use.names=TRUE,fill=TRUE)
counts_dt <- rbindlist(count_rows,use.names=TRUE,fill=TRUE)
enrich_dt <- rbindlist(enrich_all,use.names=TRUE,fill=TRUE)
enrich_dt[,supported:=is.finite(q_value)&q_value<0.05&affected_gene_count>=3L]
setorder(enrich_dt,dataset,knockout_gene,ontology,q_value,-affected_gene_count,term_id)

atomic_fwrite(affected_dt,file.path(part_dir,'phase74_affected_genes_all_v1.tsv.gz'))
atomic_fwrite(counts_dt,file.path(part_dir,'phase74_affected_gene_counts_v1.tsv'))
atomic_fwrite(enrich_dt,file.path(part_dir,'phase74_enrichment_all_v1.tsv.gz'))
atomic_fwrite(enrich_dt[supported==TRUE],file.path(part_dir,'phase74_enrichment_supported_v1.tsv'))

capture.output(sessionInfo(),file=file.path(part_dir,'R_sessionInfo_v1.txt'))
metadata <- data.table(
  status='PASS_PHASE74_ENRICHMENT_V1',input_files=length(files),
  affected_rows=nrow(affected_dt),enrichment_rows=nrow(enrich_dt),supported_rows=sum(enrich_dt$supported),
  go_annotation='org.Hs.eg.db_3.22.0',kegg_annotation='MSigDB_2026.1_Hs_C2_KEGG_LEGACY_AND_MEDICUS',
  clusterProfiler_version=as.character(packageVersion('clusterProfiler')),
  msigdbr_version=as.character(packageVersion('msigdbr')),direction_available=FALSE)
atomic_fwrite(metadata,file.path(part_dir,'phase74_run_metadata_v1.tsv'))
writeLines('PHASE74_ENRICHMENT_PASS',file.path(part_dir,'ENRICHMENT_COMPLETE'))
if (!file.rename(part_dir,output_dir)) stop('Failed to finalize output')
cat(sprintf('PASS_PHASE74_ENRICHMENT affected=%d enrichment=%d supported=%d\n',nrow(affected_dt),nrow(enrich_dt),sum(enrich_dt$supported)))
