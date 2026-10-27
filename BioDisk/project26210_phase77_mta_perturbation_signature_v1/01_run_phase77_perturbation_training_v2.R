#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(limma)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(babelgene)
  library(GenomicRanges)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args)!=2L) stop('Usage: 01_run_phase77_perturbation_training_v1.R <downloads_dir> <output_dir>')
downloads <- normalizePath(args[[1]],winslash='/',mustWork=TRUE)
out <- normalizePath(args[[2]],winslash='/',mustWork=FALSE)
part <- paste0(out,'.part')
if (dir.exists(out)||file.exists(out)||dir.exists(part)||file.exists(part)) stop('Refusing existing output')
dir.create(part,recursive=TRUE)

atomic_fwrite <- function(x,path){tmp<-paste0(path,'.part');fwrite(x,tmp,sep='\t',quote=TRUE,na='NA');if(!file.rename(tmp,path))stop('rename failed')}

map_human_ensembl <- function(ids){
  keys <- sub('\\..*$','',ids)
  unname(mapIds(org.Hs.eg.db,keys=keys,keytype='ENSEMBL',column='SYMBOL',multiVals='first'))
}

map_mouse_ensembl <- function(ids){
  keys <- sub('\\..*$','',ids)
  unname(mapIds(org.Mm.eg.db,keys=keys,keytype='ENSEMBL',column='SYMBOL',multiVals='first'))
}

map_mouse_entrez <- function(ids){
  unname(mapIds(org.Mm.eg.db,keys=as.character(ids),keytype='ENTREZID',column='SYMBOL',multiVals='first'))
}

mouse_to_human <- function(symbols){
  unique_symbols <- unique(na.omit(symbols))
  if(!length(unique_symbols)) return(data.table(symbol=character(),human_symbol=character(),support_n=integer()))
  mapped <- as.data.table(orthologs(genes=unique_symbols,species='mouse',human=FALSE))
  mapped <- mapped[!is.na(human_symbol)&human_symbol!=''&!is.na(symbol)&symbol!='']
  setorder(mapped,symbol,-support_n,human_symbol)
  mapped <- mapped[, .SD[1], by=symbol]
  mapped[,.(symbol,human_symbol=toupper(human_symbol),support_n)]
}

aggregate_counts <- function(dt,gene_col,count_cols,symbols){
  x <- copy(dt)
  x[,symbol:=toupper(symbols)]
  x <- x[!is.na(symbol)&symbol!='']
  x[,lapply(.SD,sum),by=symbol,.SDcols=count_cols]
}

run_edger <- function(count_dt,group,dataset,contrast){
  count_cols <- setdiff(names(count_dt),'symbol')
  m <- as.matrix(count_dt[,..count_cols]); storage.mode(m)<-'integer'; rownames(m)<-count_dt$symbol
  if(ncol(m)!=length(group)) stop(dataset,': group length mismatch')
  group <- factor(group,levels=c('CONTROL','PERTURBED'))
  y <- DGEList(counts=m,group=group)
  keep <- filterByExpr(y,group=group)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y,method='TMM')
  design <- model.matrix(~0+group); colnames(design)<-levels(group)
  y <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(y,design,robust=TRUE)
  qlf <- glmQLFTest(fit,contrast=c(-1,1))
  tab <- as.data.table(topTags(qlf,n=Inf,sort.by='none')$table,keep.rownames='human_symbol')
  tab[,`:=`(dataset=dataset,contrast=contrast,n_control=sum(group=='CONTROL'),n_perturbed=sum(group=='PERTURBED'))]
  setcolorder(tab,c('dataset','contrast','human_symbol','n_control','n_perturbed','logFC','logCPM','F','PValue','FDR'))
  tab[]
}

# GSE249082: direct human MTA exposure.
p249 <- file.path(downloads,'GSE249082_AC-59_RawCounts.txt.gz')
d249 <- fread(p249,skip='Geneid')
c249 <- grep('^[CM][0-9]+_',names(d249),value=TRUE)
if(length(c249)!=8L) stop('GSE249082 expected 8 count columns')
s249 <- map_human_ensembl(d249$Geneid)
a249 <- aggregate_counts(d249,'Geneid',c249,s249)
g249 <- ifelse(grepl('^C',c249),'CONTROL','PERTURBED')
de249 <- run_edger(a249,g249,'GSE249082','MTA_vs_control_human_fibroblast')

# Preserve hg38 gene coordinates for promoter ATAC audit.
split_num <- function(x,fun){v<-suppressWarnings(as.integer(strsplit(as.character(x),';',fixed=TRUE)[[1]]));if(all(is.na(v)))NA_integer_ else fun(v,na.rm=TRUE)}
coord <- data.table(human_symbol=toupper(s249),chr=sub(';.*$','',d249$Chr),strand=sub(';.*$','',d249$Strand),
                    gene_start=vapply(d249$Start,split_num,integer(1),fun=min),gene_end=vapply(d249$End,split_num,integer(1),fun=max))
coord <- coord[!is.na(human_symbol)&human_symbol!=''&grepl('^chr([0-9]+|X|Y)$',chr)&strand%chin%c('+','-')]
coord[,tss:=ifelse(strand=='+',gene_start,gene_end)]
setorder(coord,human_symbol,chr,tss); coord<-coord[, .SD[1],by=human_symbol]

# GSE276571: mouse tumor-cell Mtap knockout.
p276 <- file.path(downloads,'GSE276571_MC38_raw_counts.txt.gz')
d276 <- fread(p276)
gene276 <- if('gene'%chin%names(d276)) d276$gene else d276[[2L]]
c276 <- grep('MC38gp100_(n_DMSO|MtapKO_DMSO)_[123]$',names(d276),value=TRUE)
if(length(c276)!=6L) stop('GSE276571 expected 6 WT/KO DMSO columns')
mouse276 <- map_mouse_ensembl(gene276)
orth276 <- mouse_to_human(mouse276)
map276 <- data.table(row_id=seq_len(nrow(d276)),symbol=mouse276)[orth276,on='symbol',nomatch=0]
x276 <- cbind(data.table(row_id=seq_len(nrow(d276))),d276[,..c276])[map276,on='row_id',nomatch=0]
a276 <- x276[,lapply(.SD,sum),by=.(symbol=human_symbol),.SDcols=c276]
g276 <- ifelse(grepl('_n_DMSO_',c276),'CONTROL','PERTURBED')
de276 <- run_edger(a276,g276,'GSE276571','MtapKO_vs_WT_mouse_MC38')

# GSE298469: mouse renal IRI with MTAP inhibitor.
p298 <- file.path(downloads,'GSE298469_RAW.tar')
extract298 <- file.path(part,'GSE298469_extract')
dir.create(extract298)
untar(p298,exdir=extract298)
f298 <- list.files(extract298,pattern='GSM90150(2[0-9]|30).*Genes_ReadCount.*\\.gz$',full.names=TRUE)
get_gsm <- function(x) sub('^(GSM[0-9]+).*','\\1',basename(x))
gsm298 <- vapply(f298,get_gsm,character(1))
keep298 <- as.integer(sub('GSM','',gsm298))>=9015020 & as.integer(sub('GSM','',gsm298))<=9015030
f298<-f298[keep298];gsm298<-gsm298[keep298]
if(length(f298)!=11L) stop('GSE298469 expected 11 IRI vehicle/MTAPi files, found ',length(f298))
tables298 <- lapply(seq_along(f298),function(i){z<-fread(f298[[i]],skip='Geneid');count_col<-tail(names(z),1);z[,.(Geneid=as.character(Geneid),value=as.numeric(get(count_col)))][,sample:=gsm298[[i]]][]})
wide298 <- dcast(rbindlist(tables298),Geneid~sample,value.var='value')
c298 <- setdiff(names(wide298),'Geneid')
mouse298 <- as.character(wide298$Geneid)
orth298 <- mouse_to_human(mouse298)
map298 <- data.table(row_id=seq_len(nrow(wide298)),symbol=mouse298)[orth298,on='symbol',nomatch=0]
x298 <- cbind(data.table(row_id=seq_len(nrow(wide298))),wide298[,..c298])[map298,on='row_id',nomatch=0]
a298 <- x298[,lapply(.SD,sum),by=.(symbol=human_symbol),.SDcols=c298]
num298 <- as.integer(sub('GSM','',c298)); g298 <- ifelse(num298<=9015025,'CONTROL','PERTURBED')
de298 <- run_edger(a298,g298,'GSE298469','MTAPi_vs_vehicle_mouse_kidney_IRI')

# GSE166213: orthogonal human CD8 T-cell promoter accessibility.
p166 <- file.path(downloads,'GSE166213_RAW.tar')
extract166 <- file.path(part,'GSE166213_extract')
dir.create(extract166)
untar(p166,exdir=extract166)
f166 <- list.files(extract166,pattern='GSM50660(3[4-6]|4[0-2]).*narrowPeak.*',full.names=TRUE)
gsm166 <- vapply(f166,get_gsm,character(1))
if(length(f166)!=6L) stop('GSE166213 expected 6 mock/MTA narrowPeak files, found ',length(f166))
prom_start <- pmax(1L,coord$tss-2000L); prom_end<-coord$tss+2000L
prom <- GRanges(seqnames=coord$chr,ranges=IRanges(start=prom_start,end=prom_end),strand='*',gene=coord$human_symbol)
atac <- matrix(0,nrow=length(prom),ncol=length(f166),dimnames=list(coord$human_symbol,gsm166))
for(i in seq_along(f166)){
  z<-fread(f166[[i]],header=FALSE)
  if(ncol(z)<7L) stop('Invalid narrowPeak: ',f166[[i]])
  peaks<-GRanges(seqnames=z[[1]],ranges=IRanges(start=as.integer(z[[2]])+1L,end=as.integer(z[[3]])),signal=as.numeric(z[[7]]))
  hits<-findOverlaps(prom,peaks,ignore.strand=TRUE)
  if(length(hits)){
    sums<-tapply(mcols(peaks)$signal[subjectHits(hits)],queryHits(hits),sum)
    atac[as.integer(names(sums)),i]<-as.numeric(sums)
  }
}
keep_atac <- rowSums(atac>0)>=2L
atac_log <- log1p(atac[keep_atac,,drop=FALSE])
g166 <- factor(ifelse(as.integer(sub('GSM','',colnames(atac_log)))<=5066036,'CONTROL','PERTURBED'),levels=c('CONTROL','PERTURBED'))
design166 <- model.matrix(~0+g166);colnames(design166)<-levels(g166)
fit166 <- eBayes(contrasts.fit(lmFit(atac_log,design166),contrasts=c(-1,1)),trend=TRUE)
de166 <- as.data.table(topTable(fit166,n=Inf,sort.by='none'),keep.rownames='human_symbol')
de166[,`:=`(dataset='GSE166213',contrast='MTA_vs_mock_human_CD8_ATAC',n_control=3L,n_perturbed=3L)]
setnames(de166,c('P.Value','adj.P.Val'),c('PValue','FDR'))

# Signature construction.
setkey(de249,human_symbol);setkey(de276,human_symbol);setkey(de298,human_symbol);setkey(de166,human_symbol)
candidates <- de249[FDR<0.05 & abs(logFC)>=0.5]
candidates[,direction:=ifelse(logFC>0,'UP','DOWN')]
candidates[,weight:=pmin(abs(logFC),3)]
candidates <- merge(candidates,de276[,.(human_symbol,ko_logFC=logFC,ko_P=PValue,ko_FDR=FDR)],by='human_symbol',all.x=TRUE)
candidates <- merge(candidates,de298[,.(human_symbol,mtapi_logFC=logFC,mtapi_P=PValue,mtapi_FDR=FDR)],by='human_symbol',all.x=TRUE)
candidates <- merge(candidates,de166[,.(human_symbol,atac_logFC=logFC,atac_P=PValue,atac_FDR=FDR)],by='human_symbol',all.x=TRUE)
candidates[,`:=`(
  ko_same=!is.na(ko_P)&ko_P<0.05&sign(ko_logFC)==sign(logFC),
  mtapi_same=!is.na(mtapi_P)&mtapi_P<0.05&sign(mtapi_logFC)==sign(logFC),
  ko_opposite=!is.na(ko_P)&ko_P<0.05&sign(ko_logFC)!=sign(logFC),
  mtapi_opposite=!is.na(mtapi_P)&mtapi_P<0.05&sign(mtapi_logFC)!=sign(logFC),
  validation_testable=!is.na(ko_logFC)|!is.na(mtapi_logFC)
)]
candidates[,consensus_pass:=validation_testable&(ko_same|mtapi_same)&!(ko_opposite|mtapi_opposite)]

select_top <- function(d,n=100L){
  up <- d[direction=='UP'][order(FDR,-abs(logFC))]
  down <- d[direction=='DOWN'][order(FDR,-abs(logFC))]
  rbind(head(up,n),head(down,n),fill=TRUE)
}
direct_sig <- select_top(candidates,100L)
consensus_sig <- select_top(candidates[consensus_pass==TRUE],100L)
direct_sig[,signature_id:='DIRECT_MTA_GSE249082']
consensus_sig[,signature_id:='CONSENSUS_MTA_CROSS_CONTEXT']
signature <- rbindlist(list(direct_sig,consensus_sig),fill=TRUE)
setcolorder(signature,c('signature_id','human_symbol','direction','weight','logFC','FDR'))

counts <- signature[,.(n_genes=.N,n_up=sum(direction=='UP'),n_down=sum(direction=='DOWN')),by=signature_id]
cons_row <- counts[signature_id=='CONSENSUS_MTA_CROSS_CONTEXT']
qualified <- nrow(cons_row)==1L && cons_row$n_up>=10L && cons_row$n_down>=10L
qualification <- data.table(signature_id='CONSENSUS_MTA_CROSS_CONTEXT',n_genes=if(nrow(cons_row))cons_row$n_genes else 0L,
                            n_up=if(nrow(cons_row))cons_row$n_up else 0L,n_down=if(nrow(cons_row))cons_row$n_down else 0L,
                            qualification=if(qualified)'QUALIFIED_PRIMARY_TRANSFER' else 'CONSENSUS_SIGNATURE_NOT_QUALIFIED')

anchors <- data.table(human_symbol=c('ARG1','VEGFA','IL10','TIMP1'),expected_direction='UP')
anchors <- merge(anchors,de249[,.(human_symbol,direct_logFC=logFC,direct_P=PValue,direct_FDR=FDR)],by='human_symbol',all.x=TRUE)
anchors <- merge(anchors,de276[,.(human_symbol,ko_logFC=logFC,ko_P=PValue,ko_FDR=FDR)],by='human_symbol',all.x=TRUE)
anchors <- merge(anchors,de298[,.(human_symbol,mtapi_logFC=logFC,mtapi_P=PValue,mtapi_FDR=FDR)],by='human_symbol',all.x=TRUE)
anchors <- merge(anchors,de166[,.(human_symbol,atac_logFC=logFC,atac_P=PValue,atac_FDR=FDR)],by='human_symbol',all.x=TRUE)

atomic_fwrite(de249,file.path(part,'GSE249082_direct_MTA_edgeR_v1.tsv'))
atomic_fwrite(de276,file.path(part,'GSE276571_MtapKO_edgeR_v1.tsv'))
atomic_fwrite(de298,file.path(part,'GSE298469_MTAPi_kidney_edgeR_v1.tsv'))
atomic_fwrite(de166,file.path(part,'GSE166213_MTA_CD8_ATAC_promoter_v1.tsv'))
atomic_fwrite(candidates,file.path(part,'phase77_direct_candidates_cross_context_audit_v1.tsv'))
atomic_fwrite(signature,file.path(part,'phase77_weighted_signatures_v1.tsv'))
atomic_fwrite(counts,file.path(part,'phase77_signature_counts_v1.tsv'))
atomic_fwrite(qualification,file.path(part,'phase77_consensus_qualification_v1.tsv'))
atomic_fwrite(anchors,file.path(part,'phase77_macrophage_anchor_audit_v1.tsv'))
atomic_fwrite(coord,file.path(part,'GSE249082_hg38_gene_coordinates_v1.tsv'))

meta <- data.table(status='PASS_PHASE77_TRAINING',direct_candidates=nrow(candidates),direct_signature_genes=nrow(direct_sig),
                   consensus_signature_genes=nrow(consensus_sig),consensus_up=qualification$n_up,consensus_down=qualification$n_down,
                   consensus_qualification=qualification$qualification)
atomic_fwrite(meta,file.path(part,'run_metadata_v1.tsv'))
unlink(extract298,recursive=TRUE);unlink(extract166,recursive=TRUE)
if(!file.rename(part,out))stop('final rename failed')
cat(sprintf('{"status":"PASS_PHASE77_TRAINING","direct":%d,"consensus":%d,"qualification":"%s"}\n',nrow(direct_sig),nrow(consensus_sig),qualification$qualification))
