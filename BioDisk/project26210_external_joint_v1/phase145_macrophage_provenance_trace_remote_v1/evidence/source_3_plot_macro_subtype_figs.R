#!/usr/bin/env Rscript
# =============================================================
# 亚群级 TYMP+ <-> TAM 精准出图（读已存 cellchat_3node.rds）
#   Fig-A 三节点丰度+marker气泡  Fig-B 关键L-R聚焦气泡
#   Fig-C APP-CD74/APOE-TREM2 跨群特异性条图  Fig-D SPP1通路来源热图
# =============================================================
suppressMessages({library(CellChat); library(Seurat); library(ggplot2); library(dplyr); library(patchwork)})
PROJ<-"/mnt/BioDisk/ludingyu/project26210"
OUT<-file.path(PROJ,"results/part3_macro_subtype"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)
cc<-readRDS(file.path(OUT,"cellchat_3node.rds"))
Tp<-"TYMP-High Malignant"
TAM<-c("SPP1+ TAM","APOE+C1QC+ TAM","Other Macrophage")
df<-subsetCommunication(cc, sources.use=Tp)

## ---- Fig-B: TYMP+ -> 三类TAM 的关键L-R聚焦气泡（只留有意义的轴）----
keep_lr<-c("APP_CD74","APP_(TREM2+TYROBP)","APOE_TREM2_TYROBP","MIF_CD74_CXCR4",
           "MIF_CD74_CD44","SPP1_CD44","PROS1_AXL","C3_C3AR1","MDK_LRP1")
pB<-tryCatch(netVisual_bubble(cc, sources.use=Tp, targets.use=TAM,
      signaling=NULL, pairLR.use=NULL, remove.isolate=TRUE, thresh=0.05,
      angle.x=30, font.size=10)+ggtitle("TYMP-High -> TAM subtypes: focused L-R"),
      error=function(e){message(e);NULL})
if(!is.null(pB)) ggsave(file.path(OUT,"FigB_TYMP_to_TAM_focused.pdf"),pB,width=6.5,height=7,limitsize=FALSE)

## ---- Fig-C: 关键L-R跨所有群的特异性条图（证明 APP-CD74/APOE-TREM2 对TAM偏强）----
plot_lr_bar<-function(lr, title){
  sub<-df[df$interaction_name==lr,c("target","prob")]
  if(nrow(sub)==0) return(NULL)
  sub$target<-factor(sub$target, levels=sub$target[order(sub$prob)])
  sub$is_tam<-ifelse(grepl("SPP1\\+ TAM|APOE\\+C1QC\\+ TAM",sub$target),"target TAM",
              ifelse(grepl("Macrophage|Monocyte",sub$target),"other myeloid","non-myeloid"))
  ggplot(sub,aes(prob,target,fill=is_tam))+geom_col()+
    scale_fill_manual(values=c("target TAM"="#C0392B","other myeloid"="#E8A87C","non-myeloid"="#95A5A6"))+
    labs(x="Commun. prob.",y=NULL,fill=NULL,title=title)+theme_bw(base_size=10)+
    theme(legend.position="bottom")
}
pC1<-plot_lr_bar("APP_CD74","APP-CD74 (from TYMP-High)")
pC2<-plot_lr_bar("APOE_TREM2_TYROBP","APOE-TREM2 (from TYMP-High)")
pC<-tryCatch((pC1|pC2)+plot_layout(guides="collect")&theme(legend.position="bottom"),error=function(e)pC1)
if(!is.null(pC)) ggsave(file.path(OUT,"FigC_LR_specificity_bar.pdf"),pC,width=10,height=4.5)

## ---- Fig-D: SPP1 通路来源热图（说明 SPP1+ TAM 是SPP1主要来源，方向澄清）----
if("SPP1"%in%cc@netP$pathways){
  pdf(file.path(OUT,"FigD_SPP1_pathway_heatmap.pdf"),width=6,height=5.5)
  print(netVisual_heatmap(cc, signaling="SPP1", color.heatmap="Reds")); dev.off()
}
## APOE 通路来源热图（补充）
if("ApoE"%in%cc@netP$pathways){
  pdf(file.path(OUT,"FigD2_ApoE_pathway_heatmap.pdf"),width=6,height=5.5)
  print(netVisual_heatmap(cc, signaling="ApoE", color.heatmap="Reds")); dev.off()
}

## ---- Fig-A: 三节点 TAM marker 气泡（展示亚型定义合理）----
mac<-readRDS(file.path(PROJ,"object/02_macrophage_final_subtypes.rds"))
mac$node3<-dplyr::case_when(mac$macro_subtype=="SPP1+Macro"~"SPP1+ TAM",
                            mac$macro_subtype=="APOE+C1QC+Macro"~"APOE+C1QC+ TAM",
                            TRUE~"Other Macrophage")
Idents(mac)<-factor(mac$node3, levels=c("SPP1+ TAM","APOE+C1QC+ TAM","Other Macrophage"))
DefaultAssay(mac)<-if("RNA"%in%Assays(mac))"RNA" else DefaultAssay(mac)
mac<-NormalizeData(mac,verbose=FALSE)
mk<-c("SPP1","APOE","C1QC","C1QA","TREM2","CD74","APP","MRC1","MARCO")
mk<-intersect(mk,rownames(mac))
pA<-DotPlot(mac,features=mk)+coord_flip()+
    labs(title="TAM subtype markers")+theme_bw(base_size=10)+
    theme(axis.text.x=element_text(angle=30,hjust=1))
ggsave(file.path(OUT,"FigA_TAM_marker_dotplot.pdf"),pA,width=6,height=5)

cat("完成，输出目录:",OUT,"\n"); print(list.files(OUT,pattern="^Fig"))
