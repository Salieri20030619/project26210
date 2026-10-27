#CNV region for each sample
CRC6_CNV_regions
CRC7_CNV_regions<-read.delim("/data/huangzr/Spatial/Kras/version5/Region_devide_2/SPCNVoutputs_CRC7_clone_group_2/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat")
CRC7_CNV_regions$section<-"CRC7"


CRC6_CNV_regions<-read.delim("/data/huangzr/Spatial/Kras/version5/Region_devide_2/SPCNVoutputs_CRC6_clone_group_2/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat")
CRC6_CNV_regions$section<-"CRC6"

CRC5_CNV_regions<-read.delim("/data/huangzr/Spatial/Kras/version5/Region_devide_2/SPCNVoutputs_CRC5_clone_group_2/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat")
CRC5_CNV_regions$section<-"CRC5"

CRC4_CNV_regions<-read.delim("/data/huangzr/Spatial/Kras/version5/Region_devide_2/SPCNVoutputs_CRC4_clone_group_2/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat")
CRC4_CNV_regions$section<-"CRC4"

CRC3a_CNV_regions<-read.delim("/data/huangzr/Spatial/Kras/version5/Region_devide_2/SPCNVoutputs_CRC3a_clone_group_2/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat")
CRC3a_CNV_regions$section<-"CRC3a"

CRC3b_CNV_regions<-read.delim("/data/huangzr/Spatial/Kras/version5/Region_devide_2/SPCNVoutputs_CRC3b_clone_group_2/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat")
CRC3b_CNV_regions$section<-"CRC3b"

CRC2_CNV_regions<-read.delim("/data/huangzr/Spatial/Kras/version5/Region_devide_2/SPCNVoutputs_CRC2_clone_group_2/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat")
CRC2_CNV_regions$section<-"CRC2"

CRC1_CNV_regions<-read.delim("/data/huangzr/Spatial/Kras/version5/Region_devide_2/SPCNVoutputs_CRC1_clone_group_2/17_HMM_predHMMi6.hmm_mode-samples.pred_cnv_regions.dat")
CRC1_CNV_regions$section<-"CRC1"

Spatial_CNV_regions<-rbind(CRC1_CNV_regions,CRC2_CNV_regions,CRC3a_CNV_regions,CRC3b_CNV_regions,CRC4_CNV_regions,CRC5_CNV_regions,CRC6_CNV_regions,CRC7_CNV_regions)
table(Spatial_CNV_regions$cell_group_name)
View(Spatial_CNV_regions)
Spatial_CNV_regions$cell_group_name == "Normal_Epithelial.Normal_Epithelial_s1"
Nor_CNV_region<-paste(paste(Spatial_CNV_regions$chr[Spatial_CNV_regions$cell_group_name == "Normal_Epithelial.Normal_Epithelial_s1"],Spatial_CNV_regions$start[Spatial_CNV_regions$cell_group_name == "Normal_Epithelial.Normal_Epithelial_s1"],sep = ":"),Spatial_CNV_regions$end[Spatial_CNV_regions$cell_group_name == "Normal_Epithelial.Normal_Epithelial_s1"], sep= "-")

# move to xena browser

# CNV bins
CNV_profile<-read.csv("/data/huangzr/Spatial/Kras/version5/Fig_new/CNV/denseDataOnlyDownload.tsv.tsv",sep = "\t")
CRC5_CNV_profile<-read.csv("/data/huangzr/Spatial/Kras/version5/Fig_new/CNV/denseDataOnlyDownload_CRC5_clone.tsv",sep = "\t")
CRC5_CNV_profile$sample %in% CNV_profile$sample
CNV_profile_merge<-merge(CNV_profile,CRC5_CNV_profile,by = "sample")
nrow(CNV_profile_merge) == nrow(CNV_profile)

rownames(CNV_profile_merge)<-CNV_profile_merge$sample
CNV_profile_merge<-CNV_profile_merge[,3:ncol(CNV_profile_merge)]
CNV_profile_merge<-na.omit(CNV_profile_merge)
grep("chr6",colnames(CNV_profile_merge))
CNV_profile_merge[,c(4,26,50)]
CNV_profile_merge<-CNV_profile_merge[,c(-4,-26,-50)]
CNV_filtered_name<-colnames(CNV_profile_merge)[c(colSums(CNV_profile_merge < (-0.35) | CNV_profile_merge > (0.35))/nrow(CNV_profile_merge) > 0.05)]
CNV_profile_merge_filtered_1<-CNV_profile_merge[,colnames(CNV_profile_merge) %in% CNV_filtered_name]
colnames(CNV_profile_merge_filtered_1) <- gsub("..average.","",colnames(CNV_profile_merge_filtered_1))

Nor_CNV_region<-paste(paste(Spatial_CNV_regions$chr[Spatial_CNV_regions$cell_group_name == "Normal_Epithelial.Normal_Epithelial_s1"],Spatial_CNV_regions$start[Spatial_CNV_regions$cell_group_name == "Normal_Epithelial.Normal_Epithelial_s1"],sep = "."),Spatial_CNV_regions$end[Spatial_CNV_regions$cell_group_name == "Normal_Epithelial.Normal_Epithelial_s1"], sep= ".")
CNV_profile_merge_filtered_2 <- CNV_profile_merge_filtered_1[,!(colnames(CNV_profile_merge_filtered_1) %in% Nor_CNV_region)]

colnames(CNV_profile_merge) != "samples.y"
CNV_PCC<-cor(CNV_profile_merge[,colnames(CNV_profile_merge) != "samples.y"])
library(ggcorrplot)
rownames(CNV_PCC)<-gsub("..average.","",rownames(CNV_PCC))
colnames(CNV_PCC)<-gsub("..average.","",colnames(CNV_PCC))

library(GenomicRanges)

# 假设你的数据框名为CNV_PCC，其中列名为chromosome:start-end
# 首先，我们需要将这些列名转换为GRanges对象

rownames(CNV_PCC) <- sub("\\.", ":", rownames(CNV_PCC))

# 使用gsub函数替换剩下的.
rownames(CNV_PCC) <- gsub("\\.", "-", rownames(CNV_PCC))


colnames(CNV_PCC) <- sub("\\.", ":", colnames(CNV_PCC))

# 使用gsub函数替换剩下的.
colnames(CNV_PCC) <- gsub("\\.", "-", colnames(CNV_PCC))
colnames(CNV_PCC)
rownames(CNV_PCC) <- sub("([^-]+-[^-]+)-.*", "\\1", colnames(CNV_PCC))
CNV_PCC<-CNV_PCC[!(duplicated(sub("([^-]+-[^-]+)-.*", "\\1", colnames(CNV_PCC)))),!(duplicated(sub("([^-]+-[^-]+)-.*", "\\1", colnames(CNV_PCC))))]


saveRDS(CNV_PCC,file = "/data/huangzr/Spatial/Kras/version5/Fig_new/CNV/PCC_matrix.RDS")
h<-ComplexHeatmap::Heatmap(CNV_PCC,row_split = 10,column_split = 10,col = colorRampPalette(c("#3384BB","#A6CEE3","white","white","white",
                                                                                          "white","white","white","white","white","white","white","white","white",
                                                                                          "white","white","white","white","#FDBF6F","#FE7F02"))(100),border = TRUE)
ht <- ComplexHeatmap::draw(h)
ComplexHeatmap::row_order(ht)
pheatmap::pheatmap(CNV_PCC,scale = "none",cluster_rows = T,cluster_cols = T,breaks = seq(-1, 1, length.out = 100),
                   color = colorRampPalette(c("#3384BB","#A6CEE3","white","white","white",
                                              "white","white","white","white","white",
                                              "white","white","white","white","#FDBF6F","#FE7F02"))(100))




CNV_merge_2<-CNV_merge
CNV_merge<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/7. Spatial_CNV_region_TCGA_CNV_validate/CNV_merged_by_PCC.rds")
CNV_merge<-CNV_merge[order(CNV_merge$chr),]
CNV_merge<-paste(paste(CNV_merge$chr,CNV_merge$Start,sep = ":"),CNV_merge$End,sep = "-")
CNV_merge
CNV_merge_2
CNV_merge_frame<-data.frame(old=CNV_merge_2,new=CNV_merge)
CNV_merge_frame[!(as.character(CNV_merge_frame$old)== as.character(CNV_merge_frame$new)),]
CNV_module<-read.csv("/data/huangzr/Spatial/Kras/version5/Fig_new/denseDataOnlyDownload_update.tsv",sep = "\t")


rownames(CNV_module)<-CNV_module$sample
CNV_module<-CNV_module[,3:ncol(CNV_module)]
CNV_module<-na.omit(CNV_module)
#grep("chr6",colnames(CNV_module))
#CNV_module[,c(4,26,50)]
#CNV_module<-CNV_module[,c(-4,-26,-50)]
CNV_filtered_name<-colnames(CNV_module)[c(colSums(CNV_module < (-0.35) | CNV_module > (0.35))/nrow(CNV_module) > 0.05)]
CNV_module_filtered_1<-CNV_module[,colnames(CNV_module) %in% CNV_filtered_name]
colnames(CNV_module_filtered_1) <- gsub("..average.","",colnames(CNV_module_filtered_1))
ncol(CNV_module_filtered_1)


colnames(CNV_module_filtered_1) <- sub("\\.", ":", colnames(CNV_module_filtered_1))

# 使用gsub函数替换剩下的.
colnames(CNV_module_filtered_1) <- gsub("\\.", "-", colnames(CNV_module_filtered_1))
colnames(CNV_module_filtered_1)
CNV_module_info<-data.frame(Region=colnames(CNV_module_filtered_1),names=c("chr1q","chr1p_1","chr1p_2","chr12q_1","chr12p",
                                                          "chr12q_2","chr13q","chr14q","chr15q","chr17q",
                                                          "chr17p","chr18p","chr18q","chr20q","chr20p",
                                                          "chr22q","chr3q","chr4p","chr4q_2","chr4q_1",
                                                          "chr5p","chr7q","chr9q","chr7p","chr8p","chr8q",
                                                          "chr8p_1"))
saveRDS(CNV_module_info,file = "/data/huangzr/Spatial/Kras/version5/Fig_new/CNV_module_info.RDS")

#TCGA CNV input
CNV_module_info$Region == colnames(CNV_module_filtered_1)[1:27]
CNV_module_filtered_1$sample<-rownames(CNV_module_filtered_1)
copy_number_CNV_module<-2*(2^(CNV_module_filtered_1[,1:27]))
colnames(copy_number_CNV_module)<-paste(colnames(copy_number_CNV_module),"_adj",sep = "")
copy_number_CNV_module$sample <- rownames(copy_number_CNV_module)
colnames(CNV_module_filtered_1)[1:27]<-as.character(CNV_module_info$names)
TCGA_COAD_origin[1:5,1:5]

# 假设copy_number_CNV_module是你的数据矩阵
df <- as.data.frame(copy_number_CNV_module)  # 将矩阵转换为数据框

# 对前27列进行操作
for (i in 1:27) {
  col_name <- colnames(df)[i]
  
  # 创建新的列，将值大于等于2的放在前一列，并在列名中加后缀_gain
  df <- df %>%
    mutate(!!paste0(col_name, "_gain") := if_else(!!sym(col_name) >= 2, !!sym(col_name), NA_real_))
  
  # 创建新的列，将值小于2的放在后一列，并在列名中加后缀_loss
  df <- df %>%
    mutate(!!paste0(col_name, "_loss") := if_else(!!sym(col_name) < 2, !!sym(col_name), NA_real_))
}


# 查看结果
head(df)
df$sample<-rownames(df)
copy_number_CNV_module_sep<-df


#TCGA celltype score calculating
CNV_module_filtered_1
TCGA_COAD_origin<-readRDS("/data/huangzr/Spatial/Kras/version4/public_CRC_bulk/TCGA/TCGA_COAD_total_with_clinical.RDS")
sc_marker<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/PCC_in_real/w.ent_cluster_markers.rds")
sc_marker<-readRDS("/data/xy/Spatial_transcriptome/eWEIDE/20240326_public_Yeyouqiong/other_methods/PCC_in_real/cos_cluster_markers.rds")
View(sc_marker)
colnames(sc_marker)
sc_marker_list<-sc_marker[,c(6,7)]

library(dplyr)
# 对每一个sc_marker$cluster的factor，取出按sc_marker$p_val_adj排序的前50个基因
top_genes <- sc_marker %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 10)

clusters <- colnames(sc_marker$names)
cluster = colnames(sc_marker$names)[1]
# 为每一个cluster类别计算行均值
for (cluster in clusters) {
  TCGA_COAD_origin[[paste0(cluster, "_signature")]] <- rowMeans(
    TCGA_COAD_origin[, colnames(TCGA_COAD_origin) %in% sc_marker$names[1:50,colnames(sc_marker$names) == cluster]]
  )
}
TCGA_COAD_origin$`CD8 Trm_signature`


# CNV seperate into gain or loss group
copy_number_CNV_module
ncol(copy_number_CNV_module)
# 加载必要的包
library(dplyr)
library(tidyverse)
colMeans(TCGA_COAD_origin_tumor[,c("SPP1", "MARCO", "GPNMB", "RGCC", "MMP19", "OLR1","APOE")])


#TCGA matrix
TCGA_COAD_origin_tumor <- TCGA_COAD_origin[TCGA_COAD_origin$tumor == "Tumor",]
TCGA_COAD_origin_tumor$sample <- substr(TCGA_COAD_origin_tumor$ID, 1, 16)
TCGA_COAD_origin_tumor<-merge(TCGA_COAD_origin_tumor,CNV_module_filtered_1,by="sample")
TCGA_COAD_origin_tumor<-merge(TCGA_COAD_origin_tumor,copy_number_CNV_module_sep,by="sample")

# 找到包含"_adj"和"signature"的列
adj_cols <- grep("_adj_", colnames(TCGA_COAD_origin_tumor))
signature_cols <- grep("signature", colnames(TCGA_COAD_origin_tumor))

# 创建一个空的矩阵来保存结果
cor_matrix <- matrix(0, nrow = length(adj_cols), ncol = length(signature_cols))
rownames(cor_matrix) <- colnames(TCGA_COAD_origin_tumor)[adj_cols]
colnames(cor_matrix) <- colnames(TCGA_COAD_origin_tumor)[signature_cols]

# 对每一个包含"_adj"的列，计算其和每一个包含"signature"的列之间的皮尔森相关性
for (i in 1:length(adj_cols)) {
  for (j in 1:length(signature_cols)) {
    cor_matrix[i, j] <- cor(TCGA_COAD_origin_tumor[!(is.na(TCGA_COAD_origin_tumor[, adj_cols[i]])), adj_cols[i]], TCGA_COAD_origin_tumor[!(is.na(TCGA_COAD_origin_tumor[, adj_cols[i]])), signature_cols[j]], method = "pearson")
  }
}

for (i in 1:length(adj_cols)) {
  for (j in 1:length(signature_cols)) {
    cor_matrix[i, j] <- cor.test(TCGA_COAD_origin_tumor[!(is.na(TCGA_COAD_origin_tumor[, adj_cols[i]])), adj_cols[i]], TCGA_COAD_origin_tumor[!(is.na(TCGA_COAD_origin_tumor[, adj_cols[i]])), signature_cols[j]], method = "pearson")$p.value
  }
}
p_matrix_CRC4<-cor_matrix[c("chr20q_gain","chr20p_gain","chr7p_gain","chr7q_gain","chr17p_loss","chr17q_loss","chr8q_gain","chr8p_gain","chr1p_2_loss"),4:ncol(cor_matrix)]<0.001
cor_matrix[c("chr20q_gain","chr20p_gain","chr7p_gain","chr7q_gain","chr17p_loss","chr17q_loss",
             "chr8q_gain","chr8p_gain","chr1p_2_loss"),]>0.25
cor.test(TCGA_COAD_origin_tumor[!(is.na(TCGA_COAD_origin_tumor[, adj_cols[i]])), adj_cols[i]], TCGA_COAD_origin_tumor[!(is.na(TCGA_COAD_origin_tumor[, adj_cols[i]])), signature_cols[j]], method = "pearson")$p.value
# 打印相关性矩阵
print(cor_matrix)
cor_matrix<-cor_matrix[,4:ncol(cor_matrix)]
colnames(cor_matrix)<-gsub("_signature","",colnames(cor_matrix))
rownames(cor_matrix)<-gsub("_adj","",rownames(cor_matrix))
bk = unique(c(seq(-0.45,0.45, length=100)))
p2<-pheatmap::pheatmap(cor_matrix,color = colorRampPalette(c("#3384BB","#A6CEE3","#E4F0F6","white","white","white","white","white","white","#FEECD4","#FDBF6F","#FE7F02"))(100),cluster_rows = F,
                   cluster_cols = F,breaks = bk,border_color = "#262626",angle_col = 90,cellwidth = 15,cellheight = 15,
                   gaps_col = c(10,20,25,27,29))
p2


#subset CD4 related CNV
#View(CRC4_CNV_regions)
#c("chr14q","chr8p","chr8q","chr1p_2","chr7p","chr7q","chr14q","chr17p","chr17q","chr20p","chr20q")
#cor_matrix[c("chr14q","chr8p","chr8q","chr1p_2","chr7p","chr7q","chr14q","chr17p","chr17q","chr20p","chr20q"),]
bk = unique(c(seq(-0.4,0.4, length=100)))
display_df<-round(cor_matrix,digits = 2)
display_df[abs(display_df)<0.25] = 0
p3<-pheatmap::pheatmap(cor_matrix[c("chr20q_gain","chr20p_gain","chr7p_gain","chr7q_gain","chr17p_loss","chr17q_loss","chr8q_gain","chr8p_gain","chr1p_2_loss"),],
                       color = colorRampPalette(c("#3384BB","#A6CEE3","#E4F0F6","white","white","#FEECD4","#FDBF6F","#FE7F02"))(100),cluster_rows = F,
                       cluster_cols = F,breaks = bk,border_color = "#262626",angle_col = 90,cellwidth = 15,cellheight = 15,
                       gaps_col = c(10,20,25,27,29),gaps_row = c(4,6,8),display_numbers = display_df[c("chr20q_gain","chr20p_gain","chr7p_gain","chr7q_gain","chr17p_loss","chr17q_loss","chr8q_gain","chr8p_gain","chr1p_2_loss"),])
p3
#calculate CNV frequency for each region
CNV_module_filtered_1
CNV_Frequency<-data.frame(loss_frequency=colSums(CNV_module_filtered_1 < (-0.35))/nrow(CNV_module_filtered_1),gain_frequency=colSums(CNV_module_filtered_1 > (0.35))/nrow(CNV_module_filtered_1))
CRC4_CNV_freq<-CNV_Frequency[c("chr20q","chr7p","chr7q","chr20p","chr17p","chr17q","chr8p","chr8q","chr1p_2","chr14q"),]
CRC4_CNV_freq$region<-rownames(CRC4_CNV_freq)
CRC4_CNV_freq<-melt(CRC4_CNV_freq)
CRC4_CNV_freq$region<-factor(CRC4_CNV_freq$region,levels=rev(c("chr20q","chr20p","chr7p","chr7q","chr17p","chr17q","chr8q","chr8p","chr14q","chr1p_2")))
theme_set(theme_cowplot())
ggplot(CRC4_CNV_freq,aes(x=value,y=region,fill=variable))+geom_bar(stat = "identity",width = 0.7)+scale_fill_manual(values = c("#97C3D8","#F6B373"))


#CNV profile
CNV_clonal_region<-readRDS("/data/huangzr/Spatial/Kras/version5/Example_data/CNV_clonal_region.RDS")
clone_order<-rownames(tumor_prop_avg[c(1,6,7,9,10,2,11,8,3,4,5,12,13),c(4:6,2,3,30,29,28,19,14,10,13,21,18,31:34,20,17,1,24,11,22,9,23,26,7,12,8,15,27,16,25)])
clone_order<-gsub(" ","_",clone_order)
pheatmap::pheatmap(CNV_clonal_region[c("chr20q","chr20p","chr7p","chr7q","chr8q","chr13q","chr5p","chr8p_2","chr9q","chr12p","chr22q","chr12q_1","chr12q_2",
                                       "chr1p_2","chr8p_1"),clone_order],cluster_rows = F,cluster_cols = F,border_color = "black",
                   color = colorRampPalette(c("#5A9BC6","white","#FF942C"))(50),gaps_col = c(3),cellwidth = 10,cellheight = 10,fontsize = 8)

CNV_module_info

#Odds ratio calculating
.libPaths(c("/data/huangzr/R_package/R-3.6.1/R_package","/data/software/R/R_package","/data/xy/R/R_package/R-3.6.1/R_package","/data/software/R/R_news/R-3.6.1/lib/R/library"))

CNV_clonal_region<-CNV_clonal_region[rowSums(CNV_clonal_region == "0") != 13,]
CNV_clonal_region<-CNV_clonal_region[c("chr20q","chr20p","chr7p","chr7q","chr8q","chr13q","chr5p","chr8p_2","chr9q","chr12p","chr22q","chr12q_1","chr12q_2",
                    "chr1p_2","chr8p_1"),clone_order]

group1 <- abs(CNV_clonal_region[, 1:3])
group2 <- abs(CNV_clonal_region[, 4:ncol(CNV_clonal_region)])
ncol(group2)

OR_df<-matrix(nrow=15,ncol=4)
rownames(OR_df)<-rownames(CNV_clonal_region)
colnames(OR_df) <- c("OR","lower","upper","p")
program <- c('CNV', 'WT')
outcome <- c('Low', 'High')
library(epitools)
for (i in 0:14){
  i=i+1
  b<-sum(group1[i,])
  c<-sum(group2[i,])
  data <- matrix(c(c, b, 10-c, 3-b), nrow=2, ncol=2, byrow=TRUE)
  dimnames(data) <- list('Program'=program, 'Outcome'=outcome)
  OR_df[i,1:4]<-c(oddsratio(data,method = "small")$measure[2],oddsratio(data,method = "small")$measure[4],
                  oddsratio(data,method = "small")$measure[6],oddsratio(data,method = "small")$p.value[6])
}



OR_df[1,1]
oddsratio(data,method = "wald")
#plot
OR_df_tmp<-OR_df
OR_df<-OR_df_tmp
OR_df<-as.data.frame(OR_df)



#OR_df=OR_df[order(OR_df[,1],decreasing = F),]
OR_df$CNV<-rownames(OR_df)
OR_df$CNV <- factor(OR_df$CNV,levels = rev(OR_df$CNV))
library(ggplot2)
#OR_df$OR = OR_df$OR + 0.01
theme_set(theme_cowplot())

colnames(OR_df)[5]<-"CNV"
ggplot(OR_df, aes(x = OR, y = CNV)) + 
  geom_vline(aes(xintercept = 1), size = .25, linetype = "dashed") + 
  geom_errorbarh(aes(xmax = `upper`, xmin = `lower`), size = 0.6, 
                 height = 0.4,color="#F7A66C") +
  geom_point(size = 2,color="#F7A66C") +
  ylab("") + xlab("Odds ratio")+coord_trans(x = 'log10')+
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(0.1,1,10,100) )



# CNV status grouping 
# 定义一个函数来根据CNV score进行分组
cnv_group <- function(x) {
  ifelse(x > 0.35, "Gain", ifelse(x < -0.35, "Loss", "WT"))
}

# 对矩阵的每一行进行操作
cnv_groups <- apply(TCGA_COAD_origin_tumor[, 59436:59462], 1, function(x) {
  # 对每一行中的每一个元素使用cnv_group函数
  groups <- sapply(x, cnv_group)
  return(groups)
})

# 将结果转化为一个矩阵
cnv_groups_matrix <- matrix(unlist(cnv_groups), nrow = nrow(TCGA_COAD_origin_tumor), byrow = TRUE)

# 为新的列创建列名
colnames(TCGA_COAD_origin_tumor[, 59436:59462])
colnames(cnv_groups_matrix) <- paste0(colnames(TCGA_COAD_origin_tumor[, 59436:59462]),"_status")

# 将新的列加到原矩阵中
TCGA_COAD_origin_tumor <- cbind(TCGA_COAD_origin_tumor, cnv_groups_matrix)
colnames(TCGA_COAD_origin_tumor)[grep("status",colnames(TCGA_COAD_origin_tumor))]
grep("signature",colnames(TCGA_COAD_origin_tumor))
library(ggplot2)
library(reshape2)
library(tidyr)

TCGA_COAD_origin_tumor_2 <- TCGA_COAD_origin_tumor %>%
  mutate_at(vars(59490:59516), list(CNV_ID = ~if_else(. != "WT", paste0(colnames(TCGA_COAD_origin_tumor)[which(. != "WT") + 59489], "_", .), NA))) %>%
  pivot_longer(cols = 59490:59516, names_to = "CNV_Region", values_to = "CNV_Status", values_drop_na = TRUE)

TCGA_COAD_origin_tumor_2[,c(59402:59435)]
# 将数据转换为长格式，包括CNV区域和细胞浸润得分
long_data <- melt(TCGA_COAD_origin_tumor_2, id.vars = c("CNV_Region", "CNV_Status"), measure.vars = 59402:59435, variable.name = "Cell_Type", value.name = "Infiltration_Score")

# 分组绘制boxplot
library(ggpubr)
long_data
long_data$CNV_Region<-as.factor(long_data$CNV_Region)
i=levels(long_data$CNV_Region)[1]
for (i in levels(long_data$CNV_Region)){
  ggplot(long_data[long_data$CNV_Region == i,], aes(x = CNV_Status, y = Infiltration_Score, fill = CNV_Status)) +
    geom_boxplot() +stat_compare_means(label="p.signif",ref.group = "WT")+
    facet_wrap(~Cell_Type, scales = "free_y") +
    theme_cowplot() + scale_fill_manual(values = c("#EEA776","lightgrey","#89B8CC"))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = i, x = "CNV Status", y = "celltype signature")
  ggsave(paste("/data/huangzr/Spatial/Kras/version5/Fig_new/CNV/CNV_cellsignature_group/",i,"_signature.pdf"),width = 20,height = 20)
}
p <- ggplot(long_data, aes(x = CNV_Status, y = Infiltration_Score, fill = CNV_Status)) +
  geom_boxplot() +stat_compare_means(label="p.signif",ref.group = "WT")+
  facet_wrap(~Cell_Type, scales = "free_y") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Infiltration Scores by CNV Region and Status", x = "CNV Status", y = "Infiltration Score")
p 

# CNV co-occurance

2^0.65
colnames(copy_number_CNV_module_sep)
CNV_sep<-copy_number_CNV_module_sep[,29:82]
CNV_sep_01_matrix<-(CNV_sep > 2.549121 | CNV_sep <1.569168)
CNV_sep_01_matrix[is.na(CNV_sep_01_matrix)]<-FALSE
CNV_sep_filtered<-CNV_sep[,colSums(CNV_sep_01_matrix)/466 > 0.05]
CNV_cor_matrix<-cor(CNV_sep,use = "pairwise.complete.obs")
colnames(CNV_cor_matrix)<-gsub("_adj","",colnames(CNV_cor_matrix))
rownames(CNV_cor_matrix)<-gsub("_adj","",rownames(CNV_cor_matrix))

CNV_non_sep<-copy_number_CNV_module_sep[,1:27]
CNV_cor_matrix<-cor(CNV_non_sep)
colnames(CNV_cor_matrix)<-gsub("_adj","",colnames(CNV_cor_matrix))
rownames(CNV_cor_matrix)<-gsub("_adj","",rownames(CNV_cor_matrix))

df<-CNV_cor_matrix
cor_mat <- cor(df, use = "pairwise.complete.obs")
pval_mat <- matrix(NA, nrow(df), ncol(df))
for(i in 1:nrow(df)) {
  for(j in 1:ncol(df)) {
    temp <- cor.test(df[,i], df[,j])
    pval_mat[i,j] <- temp$p.value
  }
}

# 创建显著性星号矩阵
sig_mat <- ifelse(pval_mat < 0.001, "***", ifelse(pval_mat < 0.01, "**", ifelse(pval_mat < 0.05, "*", "")))

rownames(sig_mat) = rownames(cor_mat)
colnames(sig_mat) = colnames(cor_mat)
sig_mat<-sig_mat[c("chr20q","chr20p","chr7p","chr7q","chr17p","chr1p_2"),c("chr20q","chr20p","chr7p","chr7q","chr17p","chr1p_2")]
# 创建热图
library(pheatmap)
pheatmap(cor_mat[c("chr20q","chr20p","chr7p","chr7q","chr17p","chr1p_2"),
                                c("chr20q","chr20p","chr7p","chr7q","chr17p","chr1p_2")],
                    color = colorRampPalette(c("#3384BB","#C3DAEB","white","white","white","white","#FED9B4","#FE7F02"))(100),
                    border_color = "black",display_numbers = sig_mat)
c("#3384BB","#C3DAEB","white","white","white","white","white","white","#FED9B4","#FE7F02")
  # 添加显著性星号
heatmap + geom_text(data = melt(sig_mat), aes(Var2, Var1, label = value), size = 3)


bk = unique(c(seq(-0.9,0.9, length=100)))
pheatmap::pheatmap(CNV_cor_matrix[c("chr20q_gain","chr20p_gain","chr7p_gain","chr7q_gain","chr17p_loss","chr1p_2_loss","chr8q_gain","chr8p_gain"),c("chr20q_gain","chr20p_gain","chr7p_gain","chr7q_gain","chr17p_loss","chr1p_2_loss","chr8q_gain","chr8p_gain")],
                   scale = "none",cluster_rows = T,cluster_cols = T,breaks = seq(-1, 1, length.out = 100),
                   color = colorRampPalette(c("#3384BB","#C3DAEB","white","#FED9B4","#FE7F02"))(100),border_color = "black")
CNV_cor_matrix[c("chr20q_gain","chr20p_gain","chr7p_gain","chr7q_gain","chr17p_loss","chr17q_loss","chr8q_gain","chr8p_gain","chr1p_2_loss"),]


#CNV status
TCGA_COAD_origin_tumor$chr20q_20p_status<-paste(TCGA_COAD_origin_tumor$chr20q_status,TCGA_COAD_origin_tumor$chr20p_status)
TCGA_COAD_origin_tumor$chr20q_7p_status<-paste(TCGA_COAD_origin_tumor$chr20q_status,TCGA_COAD_origin_tumor$chr7p_status)
TCGA_COAD_origin_tumor$chr20q_7q_status<-paste(TCGA_COAD_origin_tumor$chr20q_status,TCGA_COAD_origin_tumor$chr7q_status)
TCGA_COAD_origin_tumor$chr20q_17p_status<-paste(TCGA_COAD_origin_tumor$chr20q_status,TCGA_COAD_origin_tumor$chr17p_status)
TCGA_COAD_origin_tumor$chr7p_7q_status<-paste(TCGA_COAD_origin_tumor$chr7p_status,TCGA_COAD_origin_tumor$chr7q_status)

ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_20p_status,y=LYZ))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7p_status,y=LYZ))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7q_status,y=LYZ))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("Gain Loss","WT WT","Gain WT", 'WT Loss'),],aes(x=chr20q_17p_status,y=LYZ))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
library(scales)

# 定义一个折叠函数
fold_trans <- function() {
  trans_new("fold",
            transform = function(x) ifelse(x < 30, x, ifelse(x >= 120, x - 80, 30 + (x - 30) / 9)),
            inverse = function(x) ifelse(x < 30, x, ifelse(x >= 40, x + 80, 30 + (x - 30) * 9))
  )
}
theme_set(theme_cowplot())
TCGA_COAD_origin_tumor_3<-readRDS("/data/huangzr/Spatial/Kras/version5/Fig_new/TCGA_COAD_origin_tumor_2.RDS")
load("/data/huangzr/Spatial/Kras/version5/Fig_new/work_20240411.RData")
p1<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_20p_status,y=log2(CD8A),fill=chr20q_20p_status,color=chr20q_20p_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p2<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7p_status,y=log2(CD8A),fill=chr20q_7p_status,color=chr20q_7p_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p3<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7q_status,y=log2(CD8A),fill=chr20q_7q_status,color=chr20q_7q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)
p4<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("Gain Loss","WT WT","Gain WT", 'WT Loss'),],aes(x=chr20q_17p_status,y=log2(CD8A),fill=chr20q_17p_status,color=chr20q_17p_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT WT")+scale_fill_manual(values = c("#F19695","#F4B96D","#A2C8DB","lightgrey"))+scale_color_manual(values = c("#E60013","#EF7C1A","#1E72A9","darkgrey"))+labs(x="",y="log2(TPM)")+guides(fill=F,color=F)

p1+p2+p3+p4
TCGA_COAD_origin_tumor_2[1:5,1:5]
library(scales)

TCGA_COAD_origin_tumor$chr20q_20p_status<-paste("chr20q_",TCGA_COAD_origin_tumor$chr20q_status," chr20p_",TCGA_COAD_origin_tumor$chr20p_status,sep = "")
TCGA_COAD_origin_tumor$chr20q_7p_status<-paste("chr20q_",TCGA_COAD_origin_tumor$chr20q_status," chr7p_",TCGA_COAD_origin_tumor$chr7p_status,sep = "")
TCGA_COAD_origin_tumor$chr20q_7q_status<-paste("chr20q_",TCGA_COAD_origin_tumor$chr20q_status," chr7q_",TCGA_COAD_origin_tumor$chr7q_status,sep = "")
TCGA_COAD_origin_tumor$chr20q_17p_status<-paste("chr20q_",TCGA_COAD_origin_tumor$chr20q_status," chr17p_",TCGA_COAD_origin_tumor$chr17p_status,sep = "")
TCGA_COAD_origin_tumor$chr7p_7q_status<-paste("chr7p_",TCGA_COAD_origin_tumor$chr7p_status," chr7q_",TCGA_COAD_origin_tumor$chr7q_status,sep = "")


p1<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("chr20q_Gain chr20p_Gain","chr20q_WT chr20p_WT","chr20q_Gain chr20p_WT", 'chr20q_WT chr20p_Gain'),],aes(x=chr20q_20p_status,y=CD3G,fill=chr20q_20p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr20p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p2<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("chr20q_Gain chr7p_Gain","chr20q_WT chr7p_WT","chr20q_Gain chr7p_WT", 'chr20q_WT chr7p_Gain'),],aes(x=chr20q_7p_status,y=CD3G,fill=chr20q_7p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p3<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("chr20q_Gain chr7q_Gain","chr20q_WT chr7q_WT","chr20q_Gain chr7q_WT", 'chr20q_WT chr7q_Gain'),],aes(x=chr20q_7q_status,y=CD3G,fill=chr20q_7q_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7q_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p4<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("chr20q_Gain chr17p_Loss","chr20q_WT chr17p_WT","chr20q_Gain chr17p_WT", 'chr20q_WT chr17p_Loss'),],aes(x=chr20q_17p_status,y=CD3G,fill=chr20q_17p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr17p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)

p5<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("chr20q_Gain chr20p_Gain","chr20q_WT chr20p_WT","chr20q_Gain chr20p_WT", 'chr20q_WT chr20p_Gain'),],aes(x=chr20q_20p_status,y=CD4,fill=chr20q_20p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr20p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p6<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("chr20q_Gain chr7p_Gain","chr20q_WT chr7p_WT","chr20q_Gain chr7p_WT", 'chr20q_WT chr7p_Gain'),],aes(x=chr20q_7p_status,y=CD4,fill=chr20q_7p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p7<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("chr20q_Gain chr7q_Gain","chr20q_WT chr7q_WT","chr20q_Gain chr7q_WT", 'chr20q_WT chr7q_Gain'),],aes(x=chr20q_7q_status,y=CD4,fill=chr20q_7q_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7q_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p8<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("chr20q_Gain chr17p_Loss","chr20q_WT chr17p_WT","chr20q_Gain chr17p_WT", 'chr20q_WT chr17p_Loss'),],aes(x=chr20q_17p_status,y=CD4,fill=chr20q_17p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr17p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)

p9<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("chr20q_Gain chr20p_Gain","chr20q_WT chr20p_WT","chr20q_Gain chr20p_WT", 'chr20q_WT chr20p_Gain'),],aes(x=chr20q_20p_status,y=GZMK,fill=chr20q_20p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr20p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p10<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("chr20q_Gain chr7p_Gain","chr20q_WT chr7p_WT","chr20q_Gain chr7p_WT", 'chr20q_WT chr7p_Gain'),],aes(x=chr20q_7p_status,y=GZMK,fill=chr20q_7p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p11<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("chr20q_Gain chr7q_Gain","chr20q_WT chr7q_WT","chr20q_Gain chr7q_WT", 'chr20q_WT chr7q_Gain'),],aes(x=chr20q_7q_status,y=GZMK,fill=chr20q_7q_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7q_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p12<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("chr20q_Gain chr17p_Loss","chr20q_WT chr17p_WT","chr20q_Gain chr17p_WT", 'chr20q_WT chr17p_Loss'),],aes(x=chr20q_17p_status,y=GZMK,fill=chr20q_17p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr17p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)

p13<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("chr20q_Gain chr20p_Gain","chr20q_WT chr20p_WT","chr20q_Gain chr20p_WT", 'chr20q_WT chr20p_Gain'),],aes(x=chr20q_20p_status,y=LYZ,fill=chr20q_20p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr20p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p14<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("chr20q_Gain chr7p_Gain","chr20q_WT chr7p_WT","chr20q_Gain chr7p_WT", 'chr20q_WT chr7p_Gain'),],aes(x=chr20q_7p_status,y=LYZ,fill=chr20q_7p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p15<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("chr20q_Gain chr7q_Gain","chr20q_WT chr7q_WT","chr20q_Gain chr7q_WT", 'chr20q_WT chr7q_Gain'),],aes(x=chr20q_7q_status,y=LYZ,fill=chr20q_7q_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7q_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p16<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("chr20q_Gain chr17p_Loss","chr20q_WT chr17p_WT","chr20q_Gain chr17p_WT", 'chr20q_WT chr17p_Loss'),],aes(x=chr20q_17p_status,y=LYZ,fill=chr20q_17p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr17p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)

p17<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("chr20q_Gain chr20p_Gain","chr20q_WT chr20p_WT","chr20q_Gain chr20p_WT", 'chr20q_WT chr20p_Gain'),],aes(x=chr20q_20p_status,y=C1QC,fill=chr20q_20p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr20p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p18<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("chr20q_Gain chr7p_Gain","chr20q_WT chr7p_WT","chr20q_Gain chr7p_WT", 'chr20q_WT chr7p_Gain'),],aes(x=chr20q_7p_status,y=C1QC,fill=chr20q_7p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p19<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("chr20q_Gain chr7q_Gain","chr20q_WT chr7q_WT","chr20q_Gain chr7q_WT", 'chr20q_WT chr7q_Gain'),],aes(x=chr20q_7q_status,y=C1QC,fill=chr20q_7q_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7q_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p20<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("chr20q_Gain chr17p_Loss","chr20q_WT chr17p_WT","chr20q_Gain chr17p_WT", 'chr20q_WT chr17p_Loss'),],aes(x=chr20q_17p_status,y=C1QC,fill=chr20q_17p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr17p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)

p21<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("chr20q_Gain chr20p_Gain","chr20q_WT chr20p_WT","chr20q_Gain chr20p_WT", 'chr20q_WT chr20p_Gain'),],aes(x=chr20q_20p_status,y=IGF1,fill=chr20q_20p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr20p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p22<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("chr20q_Gain chr7p_Gain","chr20q_WT chr7p_WT","chr20q_Gain chr7p_WT", 'chr20q_WT chr7p_Gain'),],aes(x=chr20q_7p_status,y=IGF1,fill=chr20q_7p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p23<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("chr20q_Gain chr7q_Gain","chr20q_WT chr7q_WT","chr20q_Gain chr7q_WT", 'chr20q_WT chr7q_Gain'),],aes(x=chr20q_7q_status,y=IGF1,fill=chr20q_7q_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7q_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p24<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("chr20q_Gain chr17p_Loss","chr20q_WT chr17p_WT","chr20q_Gain chr17p_WT", 'chr20q_WT chr17p_Loss'),],aes(x=chr20q_17p_status,y=IGF1,fill=chr20q_17p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr17p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p1+p2+p3+p4+p5+p6+p7+p8+p9+p10+p11+p12+p17+p18+p19+p20+plot_layout(ncol=4)


ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_20p_status,y=CLEC9A))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7p_status,y=CLEC9A))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7q_status,y=CLEC9A))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("Gain Loss","WT WT","Gain WT", 'WT Loss'),],aes(x=chr20q_17p_status,y=CLEC9A))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")

ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_20p_status,y=FOXP3))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7p_status,y=FOXP3))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7q_status,y=FOXP3))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("Gain Loss","WT WT","Gain WT", 'WT Loss'),],aes(x=chr20q_17p_status,y=FOXP3))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")

ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_20p_status,y=CD1C))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7p_status,y=CD1C))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7q_status,y=CD1C))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("Gain Loss","WT WT","Gain WT", 'WT Loss'),],aes(x=chr20q_17p_status,y=CD1C))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")


ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_20p_status,y=SPP1))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7p_status,y=SPP1))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7q_status,y=SPP1))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("Gain Loss","WT WT","Gain WT", 'WT Loss'),],aes(x=chr20q_17p_status,y=SPP1))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")


ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_20p_status,y=CD3G))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7p_status,y=CD3G))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("Gain Gain","WT WT","Gain WT", 'WT Gain'),],aes(x=chr20q_7q_status,y=CD3G))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("Gain Loss","WT WT","Gain WT", 'WT Loss'),],aes(x=chr20q_17p_status,y=CD3G))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "WT WT")


TCGA_COAD_origin_tumor$CNV_group<-"Remain"

TCGA_COAD_origin_tumor$other_status<-paste(paste(TCGA_COAD_origin_tumor$chr20p_status,TCGA_COAD_origin_tumor$chr7p_status),TCGA_COAD_origin_tumor$chr7q_status,TCGA_COAD_origin_tumor$chr17p_status)
TCGA_COAD_origin_tumor$CNV_group[TCGA_COAD_origin_tumor$chr20q_status == "Gain" & TCGA_COAD_origin_tumor$other_status == "WT WT WT WT"]<-"chr20q_gain_other_WT"
TCGA_COAD_origin_tumor$CNV_group[TCGA_COAD_origin_tumor$chr20q_status == "WT" & TCGA_COAD_origin_tumor$other_status == "WT WT WT WT"]<-"all_WT"
TCGA_COAD_origin_tumor$CNV_group[TCGA_COAD_origin_tumor$chr20q_status == "WT" & TCGA_COAD_origin_tumor$other_status != "WT WT WT WT"]<-"chr20q_WT_other_CNV"
p1<-ggplot(TCGA_COAD_origin_tumor,aes(x=CNV_group,y=CD3G))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "all_WT",method = "wilcox")
p2<-ggplot(TCGA_COAD_origin_tumor,aes(x=CNV_group,y=CD8A))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "all_WT",method = "wilcox")
p3<-ggplot(TCGA_COAD_origin_tumor,aes(x=CNV_group,y=CD4))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "all_WT",method = "wilcox")
p4<-ggplot(TCGA_COAD_origin_tumor,aes(x=CNV_group,y=C1QC))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "all_WT",method = "wilcox")
p1+p2+p3+p4

table(TCGA_COAD_origin_tumor$CNV_group)
TCGA_COAD_origin_tumor$CNV_group<-"Remain"
TCGA_COAD_origin_tumor$other_status<-paste(paste(TCGA_COAD_origin_tumor$chr7p_status,TCGA_COAD_origin_tumor$chr17p_status),TCGA_COAD_origin_tumor$chr20p_status,TCGA_COAD_origin_tumor$chr20q_status)
TCGA_COAD_origin_tumor$CNV_group[TCGA_COAD_origin_tumor$chr7q_status == "Gain" & TCGA_COAD_origin_tumor$other_status == "WT WT WT WT"]<-"chr7q_gain_other_WT"
TCGA_COAD_origin_tumor$CNV_group[TCGA_COAD_origin_tumor$chr7q_status == "WT" & TCGA_COAD_origin_tumor$other_status == "WT WT WT WT"]<-"all_WT"
TCGA_COAD_origin_tumor$CNV_group[TCGA_COAD_origin_tumor$chr7q_status == "WT" & TCGA_COAD_origin_tumor$other_status != "WT WT WT WT"]<-"chr7q_WT_other_CNV"
ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$CNV_group != "Remain",],aes(x=CNV_group,y=CD3G))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "all_WT",method = "wilcox")
p1<-ggplot(TCGA_COAD_origin_tumor,aes(x=CNV_group,y=CD3G))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "all_WT",method = "wilcox")
p2<-ggplot(TCGA_COAD_origin_tumor,aes(x=CNV_group,y=CD8A))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "all_WT",method = "wilcox")
p3<-ggplot(TCGA_COAD_origin_tumor,aes(x=CNV_group,y=CD4))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "all_WT",method = "wilcox")
p4<-ggplot(TCGA_COAD_origin_tumor,aes(x=CNV_group,y=C1QC))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "all_WT",method = "wilcox")
p1+p2+p3+p4

table(TCGA_COAD_origin_tumor$CNV_group)

p_unload(Seurat)
p_load(Seurat)
SpatialPlot(CRC4_st_tumor,group.by = "Region",cols = c("#ADD487","#9FC4D6","#FA7F6F"),stroke = NA,pt.size.factor = 1.2)

library(tidyr)
TCGA_COAD_origin_tumor <- TCGA_COAD_origin_tumor_3 %>% 
  pivot_wider(names_from = CNV_Region, values_from = CNV_Status)
dim(TCGA_COAD_origin_tumor)
TCGA_COAD_origin_tumor<-as.data.frame(TCGA_COAD_origin_tumor)
TCGA_COAD_origin_tumor


#TCGA type1 IFN
p1<-ggplot(TCGA_COAD_origin_tumor,aes(x=chr20q_status,y=IRF7,fill=chr20q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT")+theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c("#EBAC71","#93BDCF"))+guides(fill=F)+labs(x="chr20q status",y="Expression(TPM)",title = "IRF7")

p2<-ggplot(TCGA_COAD_origin_tumor,aes(x=chr20q_status,y=ISG15,fill=chr20q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT")+theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c("#EBAC71","#93BDCF"))+guides(fill=F)+labs(x="chr20q status",y="Expression(TPM)",title = "ISG15")

p3<-ggplot(TCGA_COAD_origin_tumor,aes(x=chr20q_status,y=STAT1,fill=chr20q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT")+theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c("#EBAC71","#93BDCF"))+guides(fill=F)+labs(x="chr20q status",y="Expression(TPM)",title = "STAT1")

p4<-ggplot(TCGA_COAD_origin_tumor,aes(x=chr20q_status,y=HLA_A,fill=chr20q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT")+theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c("#EBAC71","#93BDCF"))+guides(fill=F)+labs(x="chr20q status",y="Expression(TPM)",title = "HLA-A")

p5<-ggplot(TCGA_COAD_origin_tumor,aes(x=chr20q_status,y=B2M,fill=chr20q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT")+theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c("#EBAC71","#93BDCF"))+guides(fill=F)+labs(x="chr20q status",y="Expression(TPM)",title = "B2M")

p6<-ggplot(TCGA_COAD_origin_tumor,aes(x=chr20q_status,y=CXCL9,fill=chr20q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT")+theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c("#EBAC71","#93BDCF"))+guides(fill=F)+labs(x="chr20q status",y="Expression(TPM)",title = "CXCL9")

p7<-ggplot(TCGA_COAD_origin_tumor,aes(x=chr20q_status,y=CXCL10,fill=chr20q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT")+theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c("#EBAC71","#93BDCF"))+guides(fill=F)+labs(x="chr20q status",y="Expression(TPM)",title = "CXCL10")

p8<-ggplot(TCGA_COAD_origin_tumor,aes(x=chr20q_status,y=CXCL11,fill=chr20q_status))+geom_boxplot()+
  stat_compare_means(label="p.signif",ref.group = "WT")+theme(plot.title = element_text(hjust = 0.5))+
  scale_fill_manual(values = c("#EBAC71","#93BDCF"))+guides(fill=F)+labs(x="chr20q status",y="Expression(TPM)",title = "CXCL11")

p1+p2+p3+p4+p5+p6+p7+p8+plot_layout(ncol=4)

p1<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("chr20q_Gain chr20p_Gain","chr20q_WT chr20p_WT","chr20q_Gain chr20p_WT", 'chr20q_WT chr20p_Gain'),],aes(x=chr20q_20p_status,y=IRF7,fill=chr20q_20p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr20p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p2<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("chr20q_Gain chr7p_Gain","chr20q_WT chr7p_WT","chr20q_Gain chr7p_WT", 'chr20q_WT chr7p_Gain'),],aes(x=chr20q_7p_status,y=IRF7,fill=chr20q_7p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p3<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("chr20q_Gain chr7q_Gain","chr20q_WT chr7q_WT","chr20q_Gain chr7q_WT", 'chr20q_WT chr7q_Gain'),],aes(x=chr20q_7q_status,y=IRF7,fill=chr20q_7q_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7q_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p4<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("chr20q_Gain chr17p_Loss","chr20q_WT chr17p_WT","chr20q_Gain chr17p_WT", 'chr20q_WT chr17p_Loss'),],aes(x=chr20q_17p_status,y=IRF7,fill=chr20q_17p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr17p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)

p5<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("chr20q_Gain chr20p_Gain","chr20q_WT chr20p_WT","chr20q_Gain chr20p_WT", 'chr20q_WT chr20p_Gain'),],aes(x=chr20q_20p_status,y=STAT1,fill=chr20q_20p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr20p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p6<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("chr20q_Gain chr7p_Gain","chr20q_WT chr7p_WT","chr20q_Gain chr7p_WT", 'chr20q_WT chr7p_Gain'),],aes(x=chr20q_7p_status,y=STAT1,fill=chr20q_7p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p7<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("chr20q_Gain chr7q_Gain","chr20q_WT chr7q_WT","chr20q_Gain chr7q_WT", 'chr20q_WT chr7q_Gain'),],aes(x=chr20q_7q_status,y=STAT1,fill=chr20q_7q_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7q_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p8<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("chr20q_Gain chr17p_Loss","chr20q_WT chr17p_WT","chr20q_Gain chr17p_WT", 'chr20q_WT chr17p_Loss'),],aes(x=chr20q_17p_status,y=STAT1,fill=chr20q_17p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr17p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)

p9<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("chr20q_Gain chr20p_Gain","chr20q_WT chr20p_WT","chr20q_Gain chr20p_WT", 'chr20q_WT chr20p_Gain'),],aes(x=chr20q_20p_status,y=ISG15,fill=chr20q_20p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr20p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p10<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("chr20q_Gain chr7p_Gain","chr20q_WT chr7p_WT","chr20q_Gain chr7p_WT", 'chr20q_WT chr7p_Gain'),],aes(x=chr20q_7p_status,y=ISG15,fill=chr20q_7p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p11<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("chr20q_Gain chr7q_Gain","chr20q_WT chr7q_WT","chr20q_Gain chr7q_WT", 'chr20q_WT chr7q_Gain'),],aes(x=chr20q_7q_status,y=ISG15,fill=chr20q_7q_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7q_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p12<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("chr20q_Gain chr17p_Loss","chr20q_WT chr17p_WT","chr20q_Gain chr17p_WT", 'chr20q_WT chr17p_Loss'),],aes(x=chr20q_17p_status,y=ISG15,fill=chr20q_17p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr17p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)

p13<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_20p_status %in% c("chr20q_Gain chr20p_Gain","chr20q_WT chr20p_WT","chr20q_Gain chr20p_WT", 'chr20q_WT chr20p_Gain'),],aes(x=chr20q_20p_status,y=B2M,fill=chr20q_20p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr20p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p14<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7p_status %in% c("chr20q_Gain chr7p_Gain","chr20q_WT chr7p_WT","chr20q_Gain chr7p_WT", 'chr20q_WT chr7p_Gain'),],aes(x=chr20q_7p_status,y=B2M,fill=chr20q_7p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p15<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_7q_status %in% c("chr20q_Gain chr7q_Gain","chr20q_WT chr7q_WT","chr20q_Gain chr7q_WT", 'chr20q_WT chr7q_Gain'),],aes(x=chr20q_7q_status,y=B2M,fill=chr20q_7q_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr7q_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)
p16<-ggplot(TCGA_COAD_origin_tumor[TCGA_COAD_origin_tumor$chr20q_17p_status %in% c("chr20q_Gain chr17p_Loss","chr20q_WT chr17p_WT","chr20q_Gain chr17p_WT", 'chr20q_WT chr17p_Loss'),],aes(x=chr20q_17p_status,y=B2M,fill=chr20q_17p_status))+geom_boxplot()+stat_compare_means(label="p.signif",ref.group = "chr20q_WT chr17p_WT")+theme(axis.text.x = element_text(angle = 45,hjust = 1))+scale_fill_manual(values = c("#FB8072","#FDB462","#B3DE69","#80B1D3"))+guides(fill=F)

TCGA_COAD_origin_tumor$tumor
