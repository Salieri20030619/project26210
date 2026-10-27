#!/usr/bin/env Rscript
suppressMessages({library(Seurat); library(liana); library(dplyr); library(ggplot2); library(SingleCellExperiment)})
set.seed(42)
PROJ<-"/mnt/BioDisk/ludingyu/project26210"
OUT<-file.path(PROJ,"results/part3_liana"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)

## 1. 三节点标签
obj<-readRDS(file.path(PROJ,"object/p3_01_cellchat_input.rds"))
ct<-as.character(obj@meta.data[["cellchat_group"]]); names(ct)<-colnames(obj)
mac<-readRDS(file.path(PROJ,"object/02_macrophage_final_subtypes.rds"))
sub_map<-setNames(as.character(mac@meta.data[["macro_subtype"]]),colnames(mac))
is_mac<-ct=="Macrophage"; hit<-is_mac&(names(ct)%in%names(sub_map))
c3<-function(s)ifelse(s=="SPP1+Macro","SPP1+ TAM",ifelse(s=="APOE+C1QC+Macro","APOE+C1QC+ TAM","Other Macrophage"))
ct3<-ct; ct3[hit]<-c3(sub_map[names(ct)[hit]]); ct3[is_mac&!hit]<-"Other Macrophage"
DefaultAssay(obj)<-"RNA"

## 2. SCE
cnt<-GetAssayData(obj,assay="RNA",layer="counts"); logc<-GetAssayData(obj,assay="RNA",layer="data")
if(max(logc)==0||max(logc)>50){obj<-NormalizeData(obj,verbose=FALSE);logc<-GetAssayData(obj,layer="data")}
sce<-SingleCellExperiment(assays=list(counts=cnt,logcounts=logc))
colData(sce)$cell_type<-factor(ct3)

## 3. liana + aggregate（4 R原生方法）
res<-liana_wrap(sce, idents_col="cell_type", method=c("natmi","connectome","logfc","sca"))
saveRDS(res, file.path(OUT,"liana_res_raw.rds"))
agg<-liana_aggregate(res)
saveRDS(agg, file.path(OUT,"liana_aggregate.rds"))
write.csv(as.data.frame(lapply(agg,as.character)), file.path(OUT,"liana_aggregate_full.csv"), row.names=FALSE)
cat("agg 列:\n"); print(colnames(agg))

## 4. 稳健筛选（base R + as.character，彻底避开 %in% 报错）
d<-as.data.frame(agg)
d$source<-as.character(d$source); d$target<-as.character(d$target)
lig_col<-intersect(c("ligand.complex","ligand"),colnames(d))[1]
rec_col<-intersect(c("receptor.complex","receptor"),colnames(d))[1]
rank_col<-intersect(c("aggregate_rank","mean_rank"),colnames(d))[1]
d$lig<-as.character(d[[lig_col]]); d$rec<-as.character(d[[rec_col]])
Tp<-"TYMP-High Malignant"; TAM<-c("SPP1+ TAM","APOE+C1QC+ TAM")
foc<-d[d$source==Tp & is.element(d$target,TAM), ]
foc<-foc[order(foc[[rank_col]]), ]
cat("\n== TYMP-High -> 两类TAM 共识排名（前25，rank越小越强）==\n")
print(head(foc[,c("source","target","lig","rec",rank_col)],25))
write.csv(foc, file.path(OUT,"liana_TYMP_to_TAM.csv"), row.names=FALSE)

## 关键轴排名
cat("\n== 关键轴（APP/APOE/MIF/SPP1）在共识中的位置 ==\n")
key<-foc[grepl("APP|APOE|MIF|SPP1|C1Q|TREM2",paste(foc$lig,foc$rec)),]
print(key[,c("source","target","lig","rec",rank_col)])

## 5. 点图
p<-tryCatch(liana_dotplot(agg, source_groups=Tp, target_groups=TAM, ntop=25)+
     ggtitle("LIANA consensus (4 methods): TYMP-High -> TAM subtypes")+
     theme(axis.text.x=element_text(angle=45,hjust=1,size=9)),
     error=function(e){cat("dotplot报错:",conditionMessage(e),"\n");NULL})
if(!is.null(p)){
  ggsave(file.path(OUT,"liana_dotplot_TYMP_to_TAM.pdf"),p,width=8,height=9,limitsize=FALSE)
  ggsave(file.path(OUT,"liana_dotplot_TYMP_to_TAM.png"),p,width=8,height=9,dpi=200,limitsize=FALSE)
  cat("\n点图已存\n")
}
cat("完成 ->",OUT,"\n")
