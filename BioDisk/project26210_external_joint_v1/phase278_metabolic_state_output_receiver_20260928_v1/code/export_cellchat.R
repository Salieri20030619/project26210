p <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase278_metabolic_state_output_receiver_20260928_v1"
f <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase35_part1_c01_c04_cell_communication_v1/resource_lock_v1/CellChatDB.human_v1.rds"
db <- readRDS(f)
genes <- c("SPP1","APOE","LGALS3","TGFB1","VEGFA","IL10","IL1B","CCL2","CSF1","PDGFB","MMP9","MMP14","CTSB","CTSD","CD274")
i <- as.data.frame(db$interaction); i <- i[i$ligand %in% genes,,drop=FALSE]
write.table(i,file.path(p,"resources/CELLCHAT_OUTPUT_INTERACTIONS.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="")
co <- as.data.frame(db$complex);co$complex_name <- rownames(co)
write.table(co,file.path(p,"resources/CELLCHAT_COMPLEXES.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="")
cat("EXPORTED_INTERACTIONS",nrow(i),"\n")
