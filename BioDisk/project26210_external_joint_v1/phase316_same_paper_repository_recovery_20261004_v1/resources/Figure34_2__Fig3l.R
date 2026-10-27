# PCC between each community

total_sample_porp_immune2<-readRDS("/Path/data/total_sample_porp_immune2.RDS")
PCC_matrix <- matrix(NA, 3, 3)
for(i in 1:3) {
  for(j in 1:3) {
    temp <- cor.test(total_sample_porp_immune2[,2:4][,i], total_sample_porp_immune2[,2:4][,j],method = "pearson")
    PCC_matrix[i,j] <- temp$estimate
  }
}

rownames(PCC_matrix) = colnames(total_sample_porp_immune2[,2:4])
colnames(PCC_matrix) = colnames(total_sample_porp_immune2[,2:4])

sig_matrix <- matrix(NA, 3, 3)
for(i in 1:3) {
  for(j in 1:3) {
    temp <- cor.test(total_sample_porp_immune2[,2:4][,i], total_sample_porp_immune2[,2:4][,j],method = "pearson")
    sig_matrix[i,j] <- temp$p.value
  }
}

rownames(sig_matrix) = colnames(total_sample_porp_immune2[,2:4])
colnames(sig_matrix) = colnames(total_sample_porp_immune2[,2:4])

sig_matrix
stars_matrix <- convert_p_to_stars(sig_matrix)
stars_matrix[abs(PCC_matrix)<0.2]<-""
pheatmap::pheatmap(PCC_matrix,color = colorRampPalette(c("#3384BB","#C3DAEB","white","white","#FED9B4","#FE7F02"))(100),
                   cellwidth = 15,cellheight = 15,fontsize = 8,breaks = seq(-0.6, 0.6, length.out = 100),
                   border_color = "black",display_numbers = stars_matrix,cluster_rows = F,cluster_cols = F)
