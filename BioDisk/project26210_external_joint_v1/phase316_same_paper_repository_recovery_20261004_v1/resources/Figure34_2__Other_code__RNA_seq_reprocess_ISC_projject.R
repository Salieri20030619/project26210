# RNA-seq data reprocess
# Figure visualization
# ZuRui Huang
# Begin on May 10, 2022

# preprocess of mouse tissue RNA-seq data

rm(list=ls())
setwd("/data/huangzr/m6A_enteritis/batch1_tissue/4_final_bam")
# 1.featrueCounts
filename <- system("ls ../4_final_bam/*input.bam",intern = TRUE)
samplename <- system("ls ../4_final_bam/*final.bam|sed 's/.*\\///'|sed 's/_CleanData_rmrRNA_final.bam//'",intern = TRUE)

bamfiles <- data.frame(FileName=filename,SampleName=samplename)
sink(file = "featureCounts.log", append = FALSE, type = c("output", "message"),
     split = TRUE)
gene_count <- featureCounts(files=bamfiles$FileName,annot.ext = "/data/Genome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf",
                            isGTFAnnotationFile = TRUE,minMQS = 15, strandSpecific=2,
                            countMultiMappingReads=FALSE,isPairedEnd=TRUE,nthreads=10)
# annot.ext指定注释文件或者数据框，要求是GTF文件或者SAF格式的文件/数据框，nthreads指定线程数；
# 默认minOverlap参数值为1，allowMultiOverlap为FALSE; 默认strandSpecific是0
sink()

gene_counts <- gene_count$counts
colnames(gene_counts) <- bamfiles$SampleName
head(gene_counts)
write.csv(gene_counts,file = "gene_counts_tissue_rna_seq.csv",quote = FALSE)
gene_counts <- read.csv("gene_counts_tissue_rna_seq.csv", row.names = 1, check.names = FALSE)

sink(file = "featureCounts_rev.log", append = FALSE, type = c("output", "message"),
     split = TRUE)
gene_count_rev <- featureCounts(files=bamfiles$FileName,annot.ext = "/data/Genome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf",
                                isGTFAnnotationFile = TRUE,minMQS = 15, strandSpecific=1,
                                countMultiMappingReads=FALSE,isPairedEnd=TRUE,nthreads=10)
# annot.ext指定注释文件或者数据框，要求是GTF文件或者SAF格式的文件/数据框，nthreads指定线程数；
# 默认minOverlap参数值为1，allowMultiOverlap为FALSE; 默认strandSpecific是0
sink()

# 2. visualize the result of featureCounts
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
tmp <- read.table("featureCounts.log", sep = "\t")
assigned_ratio <- gsub(".*\\(", "", tmp[grep("Successfully assigned alignments ", tmp$V1), ])
assigned_ratio <- as.numeric(gsub("%\\).*", "", assigned_ratio))
assigned_reads <- gsub(".*: ", "", tmp[grep("Successfully assigned alignments ", tmp$V1), ])
assigned_reads <- as.numeric(gsub(" \\(.*", "", assigned_reads))
Total_alignments <- gsub(".*: ", "", tmp[grep("Total alignments : ", tmp$V1), ])
Total_alignments <- as.numeric(gsub("\\|.*", "", Total_alignments))
mapping_df <- data.frame(sample=samplename,
                         percentage=assigned_ratio[1:12], 
                         assigned_reads=assigned_reads[1:12],
                         unassigned_reads=Total_alignments[1:12]-assigned_reads[1:12])
meltdata <- reshape2::melt(mapping_df, id.vars=c("sample", "percentage"),
                           variable.name = "reads_type", value.name = "reads_number")
meltdata$reads_type <- factor(meltdata$reads_type, levels = rev(levels(meltdata$reads_type)))
setwd("/data/huangzr/m6A_enteritis/reprocess_integration")
p <- ggplot(data=meltdata) + 
  geom_bar(aes(x=sample, y=reads_number, fill=reads_type), stat="identity")+
  theme(axis.text.x  = element_text(hjust = 1, vjust = 1, angle = 30),
        legend.position = "top")+
  xlab("")+ylab("amount")+guides(fill=guide_legend(title = ""))

p + geom_point(aes(x=sample, y=percentage*2e+05))+
  geom_line(aes(x=sample, y=percentage*2e+05, group=1))+
  scale_y_continuous(sec.axis = sec_axis(~./2e+05, name="percentage",
                                         breaks=seq(0,75,25), labels=paste(seq(0,75,25), "%", sep="")))+
  geom_text(aes(x=sample, y=(percentage*2e+05)+2e+06), 
            label=paste(round(meltdata$percentage, 1), "%", sep=""))
ggsave("assigned_reads_tissue.pdf", width = 13.5, height = 7.5)
meltdata$sample <- reorder(meltdata$sample, meltdata$percentage)
p <- ggplot(data=meltdata) + 
  geom_bar(aes(x=sample, y=reads_number, fill=reads_type), stat="identity")+
  theme(axis.text.x  = element_text(hjust = 1, vjust = 1, angle = 30),
        legend.position = "top")+
  xlab("")+ylab("amount")+guides(fill=guide_legend(title = ""))

p + geom_point(aes(x=sample, y=percentage*2e+05))+
  geom_line(aes(x=sample, y=percentage*2e+05, group=1))+
  scale_y_continuous(sec.axis = sec_axis(~./2e+05, name="percentage",
                                         breaks=seq(0,75,25), labels=paste(seq(0,75,25), "%", sep="")))+
  geom_text(aes(x=sample, y=(percentage*2e+05)+2e+06), 
            label=paste(round(meltdata$percentage, 1), "%", sep=""))
ggsave("reorder_assigned_reads_tissue.pdf", width = 13.5, height = 7.5)

tmp <- read.table("featureCounts_rev.log", sep = "\t")

assigned_ratio <- gsub(".*\\(", "", tmp[grep("Successfully assigned alignments ", tmp$V1), ])
assigned_ratio <- as.numeric(gsub("%\\).*", "", assigned_ratio))
assigned_reads <- gsub(".*: ", "", tmp[grep("Successfully assigned alignments ", tmp$V1), ])
assigned_reads <- as.numeric(gsub(" \\(.*", "", assigned_reads))
Total_alignments <- gsub(".*: ", "", tmp[grep("Total alignments : ", tmp$V1), ])
Total_alignments <- as.numeric(gsub("\\|.*", "", Total_alignments))
mapping_df <- data.frame(sample=samplename,
                         percentage=assigned_ratio, 
                         assigned_reads=assigned_reads,
                         unassigned_reads=Total_alignments-assigned_reads)
meltdata <- reshape2::melt(mapping_df, id.vars=c("sample", "percentage"),
                           variable.name = "reads_type", value.name = "reads_number")
meltdata$reads_type <- factor(meltdata$reads_type, levels = rev(levels(meltdata$reads_type)))
p <- ggplot(data=meltdata) + 
  geom_bar(aes(x=sample, y=reads_number, fill=reads_type), stat="identity")+
  theme(axis.text.x  = element_text(hjust = 1, vjust = 1, angle = 30),
        legend.position = "top")+
  xlab("")+ylab("amount")+guides(fill=guide_legend(title = ""))

p + geom_point(aes(x=sample, y=percentage*2e+05))+
  geom_line(aes(x=sample, y=percentage*2e+05, group=1))+
  scale_y_continuous(sec.axis = sec_axis(~./2e+05, name="percentage",
                                         breaks=seq(0,75,25), labels=paste(seq(0,75,25), "%", sep="")))+
  geom_text(aes(x=sample, y=(percentage*2e+05)+2e+06), 
            label=paste(round(meltdata$percentage, 1), "%", sep=""))
ggsave("assigned_reads_rev_rnaseq_tissue.pdf", width = 13.5, height = 7.5)


gene_counts <- gene_counts[, grep("input", colnames(gene_counts))]
head(gene_counts)
colnames(gene_counts) <- gsub("_input", "", colnames(gene_counts))

# 3.tpm caculation
gene_length <- read.table("/data/huangzr/Genome/mm9/mm9_gene_length_max_isoform.txt",stringsAsFactors = FALSE, header = F, row.names = 1)
colnames(gene_length)<-c("longest_isoform")
gene_length <- gene_length[gene_length$longest_isoform>150, ,drop=F]
tpm_counts <- as.data.frame(gene_counts[intersect(row.names(gene_counts),row.names(gene_length)), ])
tpm_counts$gene_length <- gene_length[intersect(row.names(gene_counts),row.names(gene_length)), ]
tpm_counts <- tpm_counts/tpm_counts$gene_length
tpm_counts <- apply(tpm_counts, 2, function(x){x/sum(x)})*1000000
tpm_counts <- tpm_counts[, -ncol(tpm_counts)]
write.csv(tpm_counts,file = "tpm_counts_rna-seq_tissue.csv",quote = FALSE)
sort(rowMeans(tpm_counts),decreasing = TRUE)[5000] #37.68
sort(rowMeans(tpm_counts),decreasing = TRUE)[10000] #5.836
tpm10 <- row.names(tpm_counts[rowMeans(tpm_counts)>10,])
tpm3 <- row.names(tpm_counts[rowMeans(tpm_counts)>3,])

# 4. Gene expression distribution 
library(ggsci)
colors<-pal_locuszoom(alpha = 0.9)(7)
show_col(colors)
meltdata_TPM <- reshape2::melt(log10(tpm_counts+1), value.name="tpm")
colnames(meltdata_TPM)[1:2] <- c("gene","sample")
meltdata_TPM$sample <- gsub("_CleanData_rmrRNA$", "", meltdata_TPM$sample)
meltdata_TPM$condition <- gsub("_1$|_2$|_3$", "", meltdata_TPM$sample)
meltdata_TPM$replicate <- gsub(".*_", "", meltdata_TPM$sample)
p <- ggplot(meltdata_TPM, aes(x=tpm, color=condition))+geom_density(aes(lty=replicate)) + scale_colour_manual(values=colors[c(2,3)])
p <- p + theme(panel.grid=element_blank(),
               panel.background=element_rect(fill='transparent', color='gray'),
               axis.title = element_text(size = 22),
               axis.text = element_text(size = 20),
               legend.title = element_text(size=20, face="bold"),
               legend.text = element_text(size = 18),
               legend.key.size = unit(1.3, 'lines'))
p + xlab(expression(log[10](TPM))) 
ggsave("lgTPM_distribution_rnaseq_tissue.pdf", height = 5, width = 7.5)

meltdata_TPM <- reshape2::melt(tpm_counts, value.name="tpm")
colnames(meltdata_TPM)[1:2] <- c("gene","sample")
meltdata_TPM$condition <- gsub("_1$|_2$|_3$", "", meltdata_TPM$sample)
meltdata_TPM$replicate <- gsub(".*_", "", meltdata_TPM$sample)
p <- ggplot(meltdata_TPM, aes(x=tpm, color=condition))+geom_density(aes(lty=replicate)) + scale_colour_manual(values=colors)
p <- p + theme(panel.grid=element_blank(),
               panel.background=element_rect(fill='transparent', color='gray'),
               axis.title = element_text(size = 22),
               axis.text = element_text(size = 20),
               legend.title = element_text(size=20, face="bold"),
               legend.text = element_text(size = 18),
               legend.key.size = unit(1.3, 'lines'))
p + xlab("TPM")+xlim(c(0,20))
ggsave("TPM20_distribution.pdf", height = 5, width = 7.5)
p + xlab("TPM")+xlim(c(20,100))
ggsave("TPM100_distribution_rnaseq_tissue.pdf", height = 5, width = 7.5)

# 5. import sample information
sample_information <- readxl::read_excel("/data/huangzr/m6A_enteritis/batch1_tissue/4_final_bam/sample_information.xlsx")
coldata <- data.frame(sample_information, row.names = sample_information$`Sample name`)
colnames(gene_counts)<-gsub("_CleanData_rmrRNA","",colnames(gene_counts))
coldata <- coldata[colnames(gene_counts), ]
genematrix<- DESeqDataSetFromMatrix(countData=gene_counts,colData=coldata,design=~Condition)
# colData的行要与countData的列相对应
dim(assay(genematrix))
Filteredmatrix <- genematrix[rowMeans(counts(genematrix))>1,] # remove average count <1 gene
dim(assay(Filteredmatrix)) # 15542 genes

# 6.PCA plot
Filteredmatrix@assays@data$counts
vst <- vst(Filteredmatrix,blind=FALSE)
save(vst, file="vst_tissue_rnaseq.rdata")
load("/data/huangzr/m6A_enteritis/reprocess_integration/vst_orgnoid_rnaseq.rdata")
load("/data/huangzr/m6A_enteritis/reprocess_integration/vst_tissue_rnaseq.rdata")
# transforms the count data to the log2 scale and normalizes with respect to library size.
# blind=FALSE should be used for transforming data for downstream analysis
vsn::meanSdPlot(assay(vst),xlab="Mean counts",ylab="Standard deviation")

pcaData <- plotPCA(vst_tis,intgroup=c("Condition", "Replicate"),returnData=TRUE)
# intgroup选择colData中感兴趣的分组名,returnData选择TRUE则只返回数据不作图
percentVar <- round(100*attr(pcaData,"percentVar"))
# pcaData是数据框，其属性中记录了percentVar的值
theme_set(theme_cowplot())
pcaData$Condition<-as.character(pcaData$Condition)
pcaData$Condition[1:3] <- "Mettl3 cKO"
p <- ggplot(pcaData,aes(x=PC1,y=PC2,color=Condition,shape=as.factor(Replicate)))+geom_point(size=5)+
  xlab(paste0("PC1:",percentVar[1],"% variance"))+
  ylab(paste0("PC2:",percentVar[2],"% variance"))+
  coord_fixed(ratio = 2.3)
p +theme(panel.background=element_rect(fill='transparent', color='black',linewidth = 1.2),
         axis.line=element_line(linetype=1,color="black",size=0),
         legend.key=element_rect(fill='transparent', color='transparent'),plot.title = element_text(size=18,hjust=0.5)) + 
  labs(shape="Replicates",color="Group",title = "PCA in vivo")+
  scale_color_manual(values=c(alpha("#CC4A33",0.9),alpha("#49A5BC",0.9)))
# coord_fixed强制指定y轴和x轴坐标比值，默认为1
ggsave("PCA_invivo.pdf", height = 4, width = 5)

pcaData <- plotPCA(vst_org,intgroup=c("Condition", "Replicate"),returnData=TRUE)
# intgroup选择colData中感兴趣的分组名,returnData选择TRUE则只返回数据不作图
percentVar <- round(100*attr(pcaData,"percentVar"))
# pcaData是数据框，其属性中记录了percentVar的值
theme_set(theme_cowplot())
pcaData$Condition<-as.character(pcaData$Condition)
pcaData$Condition[1:3] <- "Mettl3 cKO"
p <- ggplot(pcaData,aes(x=PC1,y=PC2,color=Condition,shape=as.factor(Replicate)))+geom_point(size=5)+
  xlab(paste0("PC1:",percentVar[1],"% variance"))+
  ylab(paste0("PC2:",percentVar[2],"% variance"))+
  coord_fixed(ratio = 2.5)
p +theme(panel.background=element_rect(fill='transparent', color='black',linewidth = 1.2),
         axis.line=element_line(linetype=1,color="black",size=0),
         legend.key=element_rect(fill='transparent', color='transparent'),plot.title = element_text(size=18,hjust=0.5)) + 
  labs(shape="Replicates",color="Group",title = "PCA in organoid")+
  scale_color_manual(values=c(alpha("#CC4A33",0.9),alpha("#49A5BC",0.9)))
# coord_fixed强制指定y轴和x轴坐标比值，默认为1
ggsave("PCA_invivo.pdf", height = 4, width = 5)z

# 7.Differntial expressed analysis
dds <- DESeq(Filteredmatrix)
resultsNames(dds)
Mettl3_exp<-counts(dds, normalized=TRUE)["Mettl3", ]
# normalize: 根据每个基因的几何平均值得到一个伪样本 => 得到每个样本中各个基因和伪样本的比例
# 取中位数 => size factor => normalization
Mettl3_tpm<-data.frame(exp=Mettl3_exp,group=c(rep("KO",3),rep("WT",3)))
ggplot(Mettl3_tpm,aes(x=group,y=exp,fill=group))+geom_boxplot()+ scale_fill_manual(values=c("#EEA236E5","#5CB85CE5"))+
  coord_cartesian(y=c(0,430))+ylab("Normalized Counts") + labs(title = "Mettl3",x="")+
  theme(plot.title = element_text(hjust=0.5))
compare_means(exp~group,Mettl3_tpm,method = "t.test")
ggsave("Mettl3_counts_tissue.pdf", height = 5, width = 7)
norm <- counts(dds, normalized=TRUE)
head(norm)
summary(norm)
write.csv(norm, file = "norm_counts_rnaseq_tissue.csv", quote = FALSE)
upper.panel<-function(x, y){
  points(x,y, pch=".", col=c("#5BC49F"))
  r <- round(cor(x, y),digits=2)
  txt <- paste0("r = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.3, 0.9, txt)
}
read.csv("norm_counts_rnaseq_tissue.csv")
tmp <- log(norm+1, base=2)
colnames(tmp) <- gsub("", "", colnames(tmp))
pdf("pair_norm_counts.pdf", width = 10, height = 10)
pairs(tmp, lower.panel = NULL, upper.panel = upper.panel)
dev.off()

# 8.MA plot
resultsNames(dds)
result_Mettl3_KO <- results(dds, contrast=c("Condition", "KO", "WT"))
summary(result_Mettl3_KO)
write.table(cbind(row.names(result_Mettl3_KO),result_Mettl3_KO) ,"RNA_diff_tissue.txt",
            quote=FALSE,sep="\t",row.names=FALSE)
pdf("MA_plot_Mettl3_KO_tissue_rnaseq.pdf")
plotMA(result_Mettl3_KO,ylim=c(-6,6),alpha=0.01)
text(10000,4,"threshold: padj = 0.01", col = "red")
dev.off()
options(scipen = 5)
res_df_Mettl3_KO <- as.data.frame(result_Mettl3_KO)
res_df_Mettl3_KO <- na.omit(res_df_Mettl3_KO)
res_df_Mettl3_KO["Mettl3",]
res_df_Mettl3_KO$Gene <- row.names(res_df_Mettl3_KO)
res_df_Mettl3_KO$change <- ifelse(res_df_Mettl3_KO$padj<0.05 & abs(res_df_Mettl3_KO$log2FoldChange)>= 0.3,
                                  ifelse(res_df_Mettl3_KO$log2FoldChange >= 0.3,'up','down'),'not')
table(res_df_Mettl3_KO$change)
DEG_Mettl3_KO <- res_df_Mettl3_KO[res_df_Mettl3_KO$change!="not",-7]
write.csv(DEG_Mettl3_KO[order(abs(DEG_Mettl3_KO$log2FoldChange), decreasing = TRUE),],
          file="DEG_padj0.05_lfc0.3_tissue_rnaseq.csv",quote = FALSE)
DEG_Mettl3_KO_tissue<-read.csv("DEG_padj0.05_lfc0.3_tissue_rnaseq.csv",row.names = "X")
save(res_df_Mettl3_KO, file="DEG_Mettl3_KO_padj0.05_lfc0.3_tissue_rnaseq.rdata")
load("/data/huangzr/m6A_enteritis/batch1_tissue/4_final_bam/DEG_Mettl3_KO_padj0.05_lfc0.3_tissue_rnaseq.rdata")
res_df_Mettl3_KO["Tnfaip3",]
nrow(res_df_Mettl3_KO)
library(pheatmap)
annotation_row=data.frame(change=relevel(factor(DEG_Mettl3_KO_tissue$change), ref = "up"), row.names = row.names(DEG_Mettl3_KO_tissue))
pheatmap(t(scale(t(assay(vst)[row.names(DEG_Mettl3_KO_tissue),]),scale = FALSE)),
         cluster_rows =TRUE, show_colnames = TRUE, annotation_row = annotation_row,
         show_rownames = FALSE, cutree_rows = 2, cutree_cols = 2, fontsize = 15,
         breaks = seq(-2,2,length.out = 100),
         annotation_colors = list(change=c(up="#D43F3AE5", down="#46B8DAE5")),
         file="DEG_heatmap_Mettl3_KO_tissue_rnaseq.pdf", width = 6, height = 5.5, 
         color = colorRampPalette(c("#46B8DAE5","white","#D43F3AE5"))(100))
show_col(colors)
DEG_Mettl3_KO_rnaseq_tissue<-DEG_Mettl3_KO
res_df_Mettl3_KO_rnaseq_tissue<-res_df_Mettl3_KO

#orgnoid expression lfc=0.3
setwd("/data/huangzr/m6A_enteritis/batch2_orgnoid/4_final_bam")
# 1.featrueCounts
filename <- system("ls ../4_final_bam/*final_2.bam",intern = TRUE)
samplename <- system("ls ../4_final_bam/*final_2.bam|sed 's/.*\\///'|sed 's/_CleanData_rmrRNA_final_2.bam//'",intern = TRUE)

bamfiles <- data.frame(FileName=filename,SampleName=samplename)
sink(file = "featureCounts.log", append = FALSE, type = c("output", "message"),
     split = TRUE)
gene_count <- featureCounts(files=bamfiles$FileName,annot.ext = "/data/Genome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf",
                            isGTFAnnotationFile = TRUE,minMQS = 15, strandSpecific=2,
                            countMultiMappingReads=FALSE,isPairedEnd=TRUE,nthreads=10)
# annot.ext指定注释文件或者数据框，要求是GTF文件或者SAF格式的文件/数据框，nthreads指定线程数；
# 默认minOverlap参数值为1，allowMultiOverlap为FALSE; 默认strandSpecific是0
sink()

gene_counts <- gene_count$counts
colnames(gene_counts) <- bamfiles$SampleName
head(gene_counts)
write.csv(gene_counts,file = "gene_counts_orgnoid_rna_seq.csv",quote = FALSE)
gene_counts <- read.csv("gene_counts_orgnoid_rna_seq.csv", row.names = 1, check.names = FALSE)

sink(file = "featureCounts_rev.log", append = FALSE, type = c("output", "message"),
     split = TRUE)
gene_count_rev <- featureCounts(files=bamfiles$FileName,annot.ext = "/data/Genome/Mus_musculus/UCSC/mm9/Annotation/Genes/genes.gtf",
                                isGTFAnnotationFile = TRUE,minMQS = 15, strandSpecific=1,
                                countMultiMappingReads=FALSE,isPairedEnd=TRUE,nthreads=10)
# annot.ext指定注释文件或者数据框，要求是GTF文件或者SAF格式的文件/数据框，nthreads指定线程数；
# 默认minOverlap参数值为1，allowMultiOverlap为FALSE; 默认strandSpecific是0
sink()

# 2. visualize the result of featureCounts
library(ggplot2)
library(cowplot)
theme_set(theme_cowplot())
tmp <- read.table("featureCounts.log", sep = "\t")
assigned_ratio <- gsub(".*\\(", "", tmp[grep("Successfully assigned alignments ", tmp$V1), ])
assigned_ratio <- as.numeric(gsub("%\\).*", "", assigned_ratio))
assigned_reads <- gsub(".*: ", "", tmp[grep("Successfully assigned alignments ", tmp$V1), ])
assigned_reads <- as.numeric(gsub(" \\(.*", "", assigned_reads))
Total_alignments <- gsub(".*: ", "", tmp[grep("Total alignments : ", tmp$V1), ])
Total_alignments <- as.numeric(gsub("\\|.*", "", Total_alignments))
mapping_df <- data.frame(sample=samplename,
                         percentage=assigned_ratio[1:12], 
                         assigned_reads=assigned_reads[1:12],
                         unassigned_reads=Total_alignments[1:12]-assigned_reads[1:12])
meltdata <- reshape2::melt(mapping_df, id.vars=c("sample", "percentage"),
                           variable.name = "reads_type", value.name = "reads_number")
meltdata$reads_type <- factor(meltdata$reads_type, levels = rev(levels(meltdata$reads_type)))
p <- ggplot(data=meltdata) + 
  geom_bar(aes(x=sample, y=reads_number, fill=reads_type), stat="identity")+
  theme(axis.text.x  = element_text(hjust = 1, vjust = 1, angle = 30),
        legend.position = "top")+
  xlab("")+ylab("amount")+guides(fill=guide_legend(title = ""))

p + geom_point(aes(x=sample, y=percentage*2e+05))+
  geom_line(aes(x=sample, y=percentage*2e+05, group=1))+
  scale_y_continuous(sec.axis = sec_axis(~./2e+05, name="percentage",
                                         breaks=seq(0,75,25), labels=paste(seq(0,75,25), "%", sep="")))+
  geom_text(aes(x=sample, y=(percentage*2e+05)+2e+06), 
            label=paste(round(meltdata$percentage, 1), "%", sep=""))
ggsave("/data/huangzr/m6A_enteritis/reprocess_integration/assigned_reads_rnaseq_orgnoid.pdf", width = 13.5, height = 7.5)
meltdata$sample <- reorder(meltdata$sample, meltdata$percentage)
p <- ggplot(data=meltdata) + 
  geom_bar(aes(x=sample, y=reads_number, fill=reads_type), stat="identity")+
  theme(axis.text.x  = element_text(hjust = 1, vjust = 1, angle = 30),
        legend.position = "top")+
  xlab("")+ylab("amount")+guides(fill=guide_legend(title = ""))

p + geom_point(aes(x=sample, y=percentage*2e+05))+
  geom_line(aes(x=sample, y=percentage*2e+05, group=1))+
  scale_y_continuous(sec.axis = sec_axis(~./2e+05, name="percentage",
                                         breaks=seq(0,75,25), labels=paste(seq(0,75,25), "%", sep="")))+
  geom_text(aes(x=sample, y=(percentage*2e+05)+2e+06), 
            label=paste(round(meltdata$percentage, 1), "%", sep=""))
ggsave("reorder_assigned_reads_rna-seq_orgnoid.pdf", width = 13.5, height = 7.5)

tmp <- read.table("featureCounts_rev.log", sep = "\t")

assigned_ratio <- gsub(".*\\(", "", tmp[grep("Successfully assigned alignments ", tmp$V1), ])
assigned_ratio <- as.numeric(gsub("%\\).*", "", assigned_ratio))
assigned_reads <- gsub(".*: ", "", tmp[grep("Successfully assigned alignments ", tmp$V1), ])
assigned_reads <- as.numeric(gsub(" \\(.*", "", assigned_reads))
Total_alignments <- gsub(".*: ", "", tmp[grep("Total alignments : ", tmp$V1), ])
Total_alignments <- as.numeric(gsub("\\|.*", "", Total_alignments))
mapping_df <- data.frame(sample=samplename,
                         percentage=assigned_ratio, 
                         assigned_reads=assigned_reads,
                         unassigned_reads=Total_alignments-assigned_reads)
meltdata <- reshape2::melt(mapping_df, id.vars=c("sample", "percentage"),
                           variable.name = "reads_type", value.name = "reads_number")
meltdata$reads_type <- factor(meltdata$reads_type, levels = rev(levels(meltdata$reads_type)))
p <- ggplot(data=meltdata) + 
  geom_bar(aes(x=sample, y=reads_number, fill=reads_type), stat="identity")+
  theme(axis.text.x  = element_text(hjust = 1, vjust = 1, angle = 30),
        legend.position = "top")+
  xlab("")+ylab("amount")+guides(fill=guide_legend(title = ""))

p + geom_point(aes(x=sample, y=percentage*2e+05))+
  geom_line(aes(x=sample, y=percentage*2e+05, group=1))+
  scale_y_continuous(sec.axis = sec_axis(~./2e+05, name="percentage",
                                         breaks=seq(0,75,25), labels=paste(seq(0,75,25), "%", sep="")))+
  geom_text(aes(x=sample, y=(percentage*2e+05)+2e+06), 
            label=paste(round(meltdata$percentage, 1), "%", sep=""))
ggsave("assigned_reads_rev_rnaseq_orgnoid.pdf", width = 13.5, height = 7.5)

setwd("/data/huangzr/m6A_enteritis/reprocess_integration")
gene_counts <- gene_counts[, grep("input", colnames(gene_counts))]
head(gene_counts)
colnames(gene_counts) <- gsub("_input", "", colnames(gene_counts))

# 3.tpm caculation
gene_length <- read.table("/data/huangzr/Genome/mm9/mm9_gene_length_max_isoform.txt",stringsAsFactors = FALSE, header = F, row.names = 1)
colnames(gene_length)<-c("longest_isoform")
gene_length <- gene_length[gene_length$longest_isoform>150, ,drop=F]
tpm_counts <- as.data.frame(gene_counts[intersect(row.names(gene_counts),row.names(gene_length)), ])
tpm_counts$gene_length <- gene_length[intersect(row.names(gene_counts),row.names(gene_length)), ]
tpm_counts <- tpm_counts/tpm_counts$gene_length
tpm_counts <- apply(tpm_counts, 2, function(x){x/sum(x)})*1000000
tpm_counts <- tpm_counts[, -ncol(tpm_counts)]
write.csv(tpm_counts,file = "tpm_counts_orgnoid_rna-seq.csv",quote = FALSE)
sort(rowMeans(tpm_counts),decreasing = TRUE)[5000] #41.44
sort(rowMeans(tpm_counts),decreasing = TRUE)[10000] #6.881
tpm10 <- row.names(tpm_counts[rowMeans(tpm_counts)>10,])
tpm3 <- row.names(tpm_counts[rowMeans(tpm_counts)>3,])

# 4. Gene expression distribution 

meltdata_TPM_orgnoid <- reshape2::melt(log10(tpm_counts+1), value.name="tpm")
colnames(meltdata_TPM_orgnoid)[1:2] <- c("gene","sample")
meltdata_TPM_orgnoid$sample <- gsub("_CleanData_rmrRNA$", "", meltdata_TPM_orgnoid$sample)
meltdata_TPM_orgnoid$condition <- gsub("_1$|_2$|_3$", "", meltdata_TPM_orgnoid$sample)
meltdata_TPM_orgnoid$replicate <- gsub(".*_", "", meltdata_TPM_orgnoid$sample)
p <- ggplot(meltdata_TPM_orgnoid, aes(x=tpm, color=condition))+geom_density(aes(lty=replicate)) + scale_colour_manual(values=colors[c(2,3)])
p <- p + theme(panel.grid=element_blank(),
               panel.background=element_rect(fill='transparent', color='gray'),
               axis.title = element_text(size = 22),
               axis.text = element_text(size = 20),
               legend.title = element_text(size=20, face="bold"),
               legend.text = element_text(size = 18),
               legend.key.size = unit(1.3, 'lines'))
p + xlab(expression(log[10](TPM))) 
ggsave("lgTPM_distribution_rnaseq_orgnoid.pdf", height = 5, width = 7.5)

meltdata_TPM_orgnoid <- reshape2::melt(tpm_counts, value.name="tpm")
colnames(meltdata_TPM_orgnoid)[1:2] <- c("gene","sample")
meltdata_TPM_orgnoid$condition <- gsub("_1$|_2$|_3$", "", meltdata_TPM_orgnoid$sample)
meltdata_TPM_orgnoid$replicate <- gsub(".*_", "", meltdata_TPM_orgnoid$sample)
p <- ggplot(meltdata_TPM_orgnoid, aes(x=tpm, color=condition))+geom_density(aes(lty=replicate)) + scale_colour_manual(values=colors)
p <- p + theme(panel.grid=element_blank(),
               panel.background=element_rect(fill='transparent', color='gray'),
               axis.title = element_text(size = 22),
               axis.text = element_text(size = 20),
               legend.title = element_text(size=20, face="bold"),
               legend.text = element_text(size = 18),
               legend.key.size = unit(1.3, 'lines'))
p + xlab("TPM")+xlim(c(0,20))
ggsave("TPM20_distribution.pdf", height = 5, width = 7.5)
p + xlab("TPM")+xlim(c(20,100))
ggsave("TPM100_distribution_rnaseq_orgnoid.pdf", height = 5, width = 7.5)

# 5. import sample information
sample_information <- readxl::read_excel("/data/huangzr/m6A_enteritis/batch2_orgnoid/4_final_bam/sample_information.xlsx")
coldata <- data.frame(sample_information, row.names = sample_information$`Sample name`)
colnames(gene_counts)<-gsub("_CleanData_rmrRNA","",colnames(gene_counts))
colnames(gene_counts)[]
gene_counts
coldata <- coldata[colnames(gene_counts), ]
genematrix<- DESeqDataSetFromMatrix(countData=gene_counts,colData=coldata,design=~Condition)
# colData的行要与countData的列相对应
dim(assay(genematrix))
Filteredmatrix <- genematrix[rowMeans(counts(genematrix))>1,]
dim(assay(Filteredmatrix))# 15050

# 6.PCA plot
vst <- vst(Filteredmatrix,blind=FALSE)
save(vst, file="vst_orgnoid_rnaseq.rdata")
# transforms the count data to the log2 scale and normalizes with respect to library size.
# blind=FALSE should be used for transforming data for downstream analysis
vsn::meanSdPlot(assay(vst),xlab="Mean counts",ylab="Standard deviation")

pcaData <- plotPCA(vst,intgroup=c("Condition", "Replicate"),returnData=TRUE)
# intgroup选择colData中感兴趣的分组名,returnData选择TRUE则只返回数据不作图
percentVar <- round(100*attr(pcaData,"percentVar"))
# pcaData是数据框，其属性中记录了percentVar的值

p <- ggplot(pcaData,aes(x=PC1,y=PC2,color=Condition,shape=as.factor(Replicate)))+geom_point(size=3)+
  xlab(paste0("PC1:",percentVar[1],"% variance"))+
  ylab(paste0("PC2:",percentVar[2],"% variance"))+
  coord_fixed()
p + theme(panel.background=element_rect(fill='transparent', color='gray'),
          legend.key=element_rect(fill='transparent', color='transparent')) + scale_color_manual(values=c("#EEA236E5","#5CB85CE5"))
# coord_fixed强制指定y轴和x轴坐标比值，默认为1
ggsave("pca_all_sample_rna-seq_orgnoid.pdf", height = 5, width = 7)

# 7.Differntial expressed analysis
dds <- DESeq(Filteredmatrix)
resultsNames(dds)
Mettl3_exp<-counts(dds, normalized=TRUE)["Mettl3", ]
Mettl3_tpm<-data.frame(exp=Mettl3_exp,group=c(rep("KO",3),rep("WT",3)))
ggplot(Mettl3_tpm,aes(y=exp,fill=group))+geom_boxplot()+ scale_fill_manual(values=c("#EEA236E5","#5CB85CE5"))+
  coord_cartesian(y=c(0,430))+ylab("Normalized Counts") + labs(title = "Mettl3",x="")+
  theme(plot.title = element_text(hjust=0.5))
ggsave("Mettl3_counts_orgnoid.pdf", height = 5, width = 7)
norm <- counts(dds, normalized=TRUE)
head(norm)
summary(norm)
write.csv(norm, file = "norm_rnaseq_orgnoid.csv", quote = FALSE)
upper.panel<-function(x, y){
  points(x,y, pch=".", col=c("#5BC49F"))
  r <- round(cor(x, y),digits=2)
  txt <- paste0("r = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.3, 0.9, txt)
}
read.csv("norm_rnaseq_orgnoid.csv")
tmp <- log(norm+1, base=2)
colnames(tmp) <- gsub("", "", colnames(tmp))
pdf("pair_norm_counts_rnaseq_orgnoid.pdf", width = 10, height = 10)
pairs(tmp, lower.panel = NULL, upper.panel = upper.panel)
dev.off()

# 8.MA plot
resultsNames(dds)
result_Mettl3_KO <- results(dds, contrast=c("Condition", "KO", "WT"))
summary(result_Mettl3_KO)
write.table(cbind(row.names(result_Mettl3_KO),result_Mettl3_KO) ,"RNA_diff_orgnoid.txt",
            quote=FALSE,sep="\t",row.names=FALSE)
pdf("MA_plot_Mettl3_KO_orgnoid_rnaseq.pdf")
plotMA(result_Mettl3_KO,ylim=c(-6,6),alpha=0.01)
text(10000,4,"threshold: padj = 0.01", col = "red")
dev.off()
options(scipen = 5)
res_df_Mettl3_KO <- as.data.frame(result_Mettl3_KO)
res_df_Mettl3_KO <- na.omit(res_df_Mettl3_KO)
res_df_Mettl3_KO["Mettl3",]
res_df_Mettl3_KO$Gene <- row.names(res_df_Mettl3_KO)
res_df_Mettl3_KO$change <- ifelse(res_df_Mettl3_KO$pvalue<0.05 & abs(res_df_Mettl3_KO$log2FoldChange)>= 0.3,
                                  ifelse(res_df_Mettl3_KO$log2FoldChange >= 0.3,'up','down'),'not')
table(res_df_Mettl3_KO$change)
DEG_Mettl3_KO <- res_df_Mettl3_KO[res_df_Mettl3_KO$change!="not",-7]
write.csv(DEG_Mettl3_KO[order(abs(DEG_Mettl3_KO$log2FoldChange), decreasing = TRUE),],
          file="DEG_Mettl3_KO_pval0.05_lfc0.3_orgnoid_rnaseq.csv",quote = FALSE)
save(res_df_Mettl3_KO, file="DEG_Mettl3_KO_pval0.05_lfc0.3_orgnoid_rnaseq.rdata")
load("DEG_Mettl3_KO_padj0.05_lfc0.3_orgnoid_rnaseq.rdata")
res_df_Mettl3_KO
dev.new()
library(pheatmap)
annotation_row=data.frame(change=relevel(factor(DEG_Mettl3_KO_orgnoid$change), ref = "up"), row.names = row.names(DEG_Mettl3_KO_orgnoid))
pheatmap(t(scale(t(assay(vst)[row.names(DEG_Mettl3_KO_tissue),]),scale = FALSE)),
         cluster_rows =TRUE, show_colnames = TRUE, annotation_row = annotation_row,
         show_rownames = FALSE, cutree_rows = 2, cutree_cols = 2, fontsize = 15,
         breaks = seq(-2,2,length.out = 100), 
         annotation_colors = list(change=c(up="#D43F3AE5", down="#46B8DAE5")),
         width = 6, height = 5.5, 
         color = colorRampPalette(c("#46B8DAE5","white","#D43F3AE5"))(100))

DEG_Mettl3_KO_rnaseq_orgnoid<-DEG_Mettl3_KO
res_df_Mettl3_KO_rnaseq_orgnoid<-res_df_Mettl3_KO
table(rownames(DEG_Mettl3_KO_rnaseq_orgnoid) %in% rownames(DEG_Mettl3_KO_rnaseq_tissue))
res_df_Mettl3_KO_rnaseq_tissue
#integrate
res_df_Mettl3_KO_rnaseq_tissue<-res_df_Mettl3_KO_rnaseq_tissue[,c(1,2,6,8)]
res_df_Mettl3_KO_rnaseq_orgnoid<-res_df_Mettl3_KO_rnaseq_orgnoid[,c(1,2,6,8)]
colnames(res_df_Mettl3_KO_rnaseq_tissue)<-c("baseMean_tis","lfc_tis","padj_tis","change_tis")
colnames(res_df_Mettl3_KO_rnaseq_orgnoid)<-c("baseMean_org","lfc_org","padj_org","change_org")
res_df_Mettl3_KO_rnaseq_orgnoid$gene=rownames(res_df_Mettl3_KO_rnaseq_orgnoid)
res_df_Mettl3_KO_rnaseq_tissue$gene=rownames(res_df_Mettl3_KO_rnaseq_tissue)