p <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase316_same_paper_repository_recovery_20261004_v1"
stopifnot(system("/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh")==0)
x <- readRDS(file.path(p,"downloads/Example_Data_Human_CRC_ST_data.rds"))
cat("CLASS",class(x),"\n");cat("DIM",dim(x),"\n");cat("SLOTS",slotNames(x),"\n")
if (isS4(x) && "meta.data" %in% slotNames(x)) {
 m <- as.data.frame(slot(x,"meta.data"));cat("META_COLUMNS",names(m),"\n");print(utils::head(m[,seq_len(min(12,ncol(m))),drop=FALSE],3))
 write.table(cbind(barcode=rownames(m),m),file.path(p,"tables/EXAMPLE_METADATA.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
 cat("IMAGES",names(slot(x,"images")),"\n")
 for (nm in names(slot(x,"images"))) {
  im <- slot(x,"images")[[nm]];cat("IMAGE",nm,"CLASS",class(im),"SLOTS",slotNames(im),"\n")
  if("coordinates" %in% slotNames(im)) {co <- as.data.frame(slot(im,"coordinates"));print(head(co,3));write.table(cbind(barcode=rownames(co),co),file.path(p,paste0("tables/EXAMPLE_COORDINATES_",nm,".tsv")),sep="\t",quote=FALSE,row.names=FALSE)}
 }
 cat("ASSAYS",names(slot(x,"assays")),"\n")
}
writeLines(capture.output(sessionInfo()),file.path(p,"resources/R_SESSION.txt"))
