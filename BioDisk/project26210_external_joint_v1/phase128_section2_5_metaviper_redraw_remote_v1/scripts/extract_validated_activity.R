options(digits=17)
w <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase128_section2_5_metaviper_redraw_remote_v1"
base <- "/mnt/BioDisk/ludingyu/project26210/results/part1_rebuild_v1/08_metaviper_v1"
a <- readRDS(file.path(base,"rds/metaviper_activity_entrez_by_cell_v1.rds"))
m <- as.data.frame(readRDS(file.path(base,"rds/metaviper_cell_metadata_v1.rds")))
mp <- read.delim(file.path(base,"tables/metaviper_regulator_id_map_v1.tsv"),colClasses="character")
genes <- c("CASP9","CHMP1A","PRSS8","PXN","HIF1A")
ids <- mp$entrez_id[match(genes,mp$symbol)]
stopifnot(identical(dim(a),c(5566L,11921L)),identical(colnames(a),m$exact_cell_id),length(unique(m$curated_patient_id))==20L,all(ids %in% rownames(a)))
d <- data.frame(exact_cell_id=colnames(a),curated_patient_id=m$curated_patient_id,as.data.frame(t(a[ids,,drop=FALSE])),check.names=FALSE)
names(d)[3:7] <- genes
stopifnot(all(is.finite(as.matrix(d[,genes]))))
write.table(d,file.path(w,"data/SELECTED_METAVIPER_ACTIVITY.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(w,"qc/R_SESSION_INFO.txt"))
cat("EXTRACTED_11921_CELLS_5_FEATURES\n")
