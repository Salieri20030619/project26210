#!/usr/bin/env Rscript
# =============================================================
# 亚群级 CellChat 重跑：七群版 vs 三节点版 对比
#   输入: p3_01_cellchat_input.rds(含 TYMP+/TYMP- 恶性等12群) + 02_macrophage_final_subtypes.rds(macro_subtype)
#   输出: results/part3_macro_subtype/ 下两版通讯图 + 对比CSV
# =============================================================
suppressMessages({library(Seurat); library(CellChat); library(patchwork); library(dplyr)})
set.seed(42)
PROJ <- "/mnt/BioDisk/ludingyu/project26210"
IN   <- file.path(PROJ,"object/p3_01_cellchat_input.rds")
MAC  <- file.path(PROJ,"object/02_macrophage_final_subtypes.rds")
OUT  <- file.path(PROJ,"results/part3_macro_subtype"); dir.create(OUT,showWarnings=FALSE,recursive=TRUE)

## ---------- 0. 读入 + 确认标签列 ----------
obj <- readRDS(IN)
cat("== p3_01 输入对象 meta 列 ==\n"); print(colnames(obj@meta.data))
# 自动探测细胞类型列（含 Macrophage 且含 TYMP+ Malignant 的那一列）
ct_col <- NULL
for(c in colnames(obj@meta.data)){
  v <- as.character(obj@meta.data[[c]])
  if(any(grepl("Macrophage",v)) && any(grepl("TYMP",v))){ ct_col <- c; break }
}
if(is.null(ct_col)) stop("未找到同时含 Macrophage 与 TYMP 的细胞类型列，请手动指定 ct_col")
cat(">> 使用细胞类型列:", ct_col, "\n")
cat(">> 原始类别:\n"); print(table(obj@meta.data[[ct_col]]))

## ---------- 1. 把 macro_subtype 标签 merge 回来 ----------
mac <- readRDS(MAC)
sub_map <- setNames(as.character(mac@meta.data[["macro_subtype"]]), colnames(mac))
ct <- as.character(obj@meta.data[[ct_col]])
names(ct) <- colnames(obj)
is_mac <- ct=="Macrophage"
cat("\n输入对象 Macrophage 数:", sum(is_mac),
    " | 能匹配到亚型的:", sum(names(ct)[is_mac] %in% names(sub_map)), "\n")

## ---------- 2. 构建两版标签向量 ----------
# 七群版
ct7 <- ct
hit <- is_mac & (colnames(obj) %in% names(sub_map))
ct7[hit] <- sub_map[colnames(obj)[hit]]
# 三节点版：SPP1+ / APOE+C1QC+ / 其余Other
collapse3 <- function(s){
  ifelse(s=="SPP1+Macro","SPP1+ TAM",
    ifelse(s=="APOE+C1QC+Macro","APOE+C1QC+ TAM",
      ifelse(grepl("Macro$",s),"Other Macrophage", s)))
}
ct3 <- ct7; ct3[hit] <- collapse3(sub_map[colnames(obj)[hit]])
# 未匹配到亚型的 Macrophage 兜底并入 Other
ct7[is_mac & !hit] <- "Other Macrophage"
ct3[is_mac & !hit] <- "Other Macrophage"

cat("\n== 七群版类别 ==\n"); print(table(ct7))
cat("\n== 三节点版类别 ==\n"); print(table(ct3))

## ---------- 3. CellChat 运行函数 ----------
run_cc <- function(labels, tag){
  cat("\n########## 运行", tag, "##########\n")
  meta <- data.frame(labels=factor(labels), row.names=colnames(obj))
  ay <- if("RNA"%in%Assays(obj))"RNA" else DefaultAssay(obj)
  data.input <- tryCatch(GetAssayData(obj, assay=ay, layer="data"),
                         error=function(e) GetAssayData(obj, assay=ay, slot="data"))
  if(max(data.input)>50){ cat(">> 检测到非log数据，NormalizeData\n"); obj2<-NormalizeData(obj,verbose=FALSE); data.input<-GetAssayData(obj2,assay=ay,layer="data") }
  cc <- createCellChat(object=data.input, meta=meta, group.by="labels")
  cc@DB <- CellChatDB.human
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type="triMean", seed.use=42)
  cc <- filterCommunication(cc, min.cells=10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  saveRDS(cc, file.path(OUT, paste0("cellchat_",tag,".rds")))
  cc
}

cc7 <- run_cc(ct7, "7group")
cc3 <- run_cc(ct3, "3node")

## ---------- 4. 出图：聚焦 TYMP+ 与各 TAM 亚型 ----------
plot_set <- function(cc, tag, tam_levels){
  lv <- levels(cc@idents)
  tam <- intersect(tam_levels, lv)
  Tp <- "TYMP-High Malignant"
  if(!Tp %in% lv){ cat("警告:",tag,"缺 TYMP-High Malignant\n"); return(invisible()) }
  # (a) TYMP+ -> 各TAM 的配受体气泡
  p1 <- tryCatch(netVisual_bubble(cc, sources.use=Tp, targets.use=tam,
          remove.isolate=TRUE, thresh=0.05, angle.x=45, font.size=8)+
          ggtitle(paste0(tag,": TYMP-High -> TAM subtypes")),
        error=function(e){message(e);NULL})
  if(!is.null(p1)) ggsave(file.path(OUT,paste0(tag,"_A_TYMP_to_TAM.pdf")),p1,width=7,height=8,limitsize=FALSE)
  # (b) 各TAM -> TYMP+
  p2 <- tryCatch(netVisual_bubble(cc, sources.use=tam, targets.use=Tp,
          remove.isolate=TRUE, thresh=0.05, angle.x=45, font.size=8)+
          ggtitle(paste0(tag,": TAM subtypes -> TYMP-High")),
        error=function(e){message(e);NULL})
  if(!is.null(p2)) ggsave(file.path(OUT,paste0(tag,"_B_TAM_to_TYMP.pdf")),p2,width=7,height=8,limitsize=FALSE)
  # (c) SPP1/MIF/APP/APOE 通路热图
  for(pw in intersect(c("SPP1","MIF","APP","APOE"), cc@netP$pathways)){
    pdf(file.path(OUT,paste0(tag,"_path_",pw,"_heatmap.pdf")),width=6,height=5.5)
    print(netVisual_heatmap(cc, signaling=pw, color.heatmap="Reds")); dev.off()
  }
  # (d) 交互强度矩阵（导出数值供对比）
  w <- cc@net$weight
  write.csv(w, file.path(OUT,paste0(tag,"_interaction_weight_matrix.csv")))
}
plot_set(cc7,"7group", c("SPP1+Macro","APOE+C1QC+Macro","CX3CR1+Macro","IL1B_IL8+Macro","HLA-II+Macro","CCL3+Macro","Cycling+Macro","Other Macrophage"))
plot_set(cc3,"3node",  c("SPP1+ TAM","APOE+C1QC+ TAM","Other Macrophage"))

## ---------- 5. 两版对比表：TYMP+ 对每类TAM的总互作强度 ----------
extract_TYMP_to <- function(cc,tag){
  w<-cc@net$weight; Tp<-"TYMP-High Malignant"
  if(!Tp%in%rownames(w)) return(NULL)
  data.frame(version=tag, target=colnames(w), TYMP_to_target=w[Tp,], target_to_TYMP=w[,Tp])
}
cmp <- rbind(extract_TYMP_to(cc7,"7group"), extract_TYMP_to(cc3,"3node"))
write.csv(cmp, file.path(OUT,"CMP_TYMP_vs_TAM_strength.csv"), row.names=FALSE)
cat("\n===== TYMP+ 对各群互作强度（两版对比）=====\n"); print(cmp)
cat("\n全部输出 ->", OUT, "\n")
