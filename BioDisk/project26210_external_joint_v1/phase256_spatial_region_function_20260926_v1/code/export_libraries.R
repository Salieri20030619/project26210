args<-commandArgs(TRUE); p<-args[1]
suppressPackageStartupMessages({library(AnnotationDbi);library(org.Hs.eg.db);library(GO.db)})
db<-org.Hs.eg.db
anno<-select(db,keys=keys(db,keytype="ENTREZID"),columns=c("SYMBOL","GOALL","ONTOLOGYALL"),keytype="ENTREZID")
anno<-unique(anno[!is.na(anno$SYMBOL)&!is.na(anno$GOALL)&anno$ONTOLOGYALL=="BP",c("GOALL","SYMBOL")])
ids<-sort(unique(anno$GOALL)); info<-select(GO.db,keys=ids,columns="TERM",keytype="GOID")
terms<-setNames(info$TERM,info$GOID); sets<-split(anno$SYMBOL,anno$GOALL)
con<-file(file.path(p,"resources","GO_BP_local.gmt"),"wt")
for(id in sort(names(sets)))writeLines(paste(c(paste0(id,"__",gsub("[\t\r\n]"," ",terms[id])),id,sort(unique(sets[[id]]))),collapse="\t"),con)
close(con)
write.table(anno,file.path(p,"resources","GO_BP_symbol_mapping.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c(capture.output(sessionInfo()),capture.output(as.data.frame(metadata(db))),capture.output(as.data.frame(metadata(GO.db)))),file.path(p,"resources","R_ANNOTATION_VERSIONS.txt"))
cat("GO_BP_SETS",length(sets),"\n")
