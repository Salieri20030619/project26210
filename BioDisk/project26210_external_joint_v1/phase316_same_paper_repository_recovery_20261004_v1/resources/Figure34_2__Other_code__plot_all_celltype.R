library(gridExtra)
library(grid)
library(gtable)
width = 21
height = 29.7


plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(UCASpatial[[1]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("UCASpatial_CRC1_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(UCASpatial[[4]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("UCASpatial_CRC2_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(UCASpatial[[2]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("UCASpatial_CRC3a_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(UCASpatial[[6]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("UCASpatial_CRC4_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(UCASpatial[[5]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("UCASpatial_CRC5_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(UCASpatial[[7]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("UCASpatial_CRC7_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()

plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(UCASpatial[[8]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("UCASpatial_CRC7_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
library(gridExtra)
library(grid)
library(gtable)
width = 21
height = 29.7


plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(SPOTlight[[1]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("SPOTlight_CRC1_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(SPOTlight[[4]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("SPOTlight_CRC2_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(SPOTlight[[2]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("SPOTlight_CRC3a_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(SPOTlight[[6]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("SPOTlight_CRC4_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(SPOTlight[[5]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("SPOTlight_CRC5_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(SPOTlight[[7]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("SPOTlight_CRC7_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()

plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(SPOTlight[[8]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("SPOTlight_CRC7_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
library(gridExtra)
library(grid)
library(gtable)
width = 21
height = 29.7


plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(CARD[[1]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("CARD_CRC1_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(CARD[[4]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("CARD_CRC2_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(CARD[[2]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("CARD_CRC3a_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(CARD[[6]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("CARD_CRC4_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(CARD[[5]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("CARD_CRC5_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(CARD[[7]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("CARD_CRC7_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()

plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(CARD[[8]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("CARD_CRC7_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
library(gridExtra)
library(grid)
library(gtable)
width = 21
height = 29.7


plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(steroscope[[1]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("steroscope_CRC1_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(steroscope[[4]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("steroscope_CRC2_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(steroscope[[2]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("steroscope_CRC3a_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(steroscope[[6]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("steroscope_CRC4_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(steroscope[[5]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("steroscope_CRC5_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(steroscope[[7]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("steroscope_CRC7_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()

plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(steroscope[[8]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("steroscope_CRC7_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
library(gridExtra)
library(grid)
library(gtable)
width = 21
height = 29.7


plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(C2L[[1]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("C2L_CRC1_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(C2L[[4]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("C2L_CRC2_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(C2L[[2]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("C2L_CRC3a_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(C2L[[6]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("C2L_CRC4_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(C2L[[5]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("C2L_CRC5_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(C2L[[7]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("C2L_CRC7_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()

plots <- lapply(celltype_all, function(ct) {
  p <- SpatialFeaturePlot(C2L[[8]], features = c(ct), stroke = NA,min.cutoff = 0.03,alpha = c(0.1,1)) +
    theme(legend.position = "right",legend.title = element_blank()) + 
    scale_fill_viridis(direction = 1)+ggtitle(c(ct))+theme(plot.title = element_text(hjust=0.5))
  return(p)
})
pdf("C2L_CRC7_cut0.03.pdf", width = width, height = height)
grid.arrange(grobs=plots, ncol = 6,vp=viewport(width=0.85, height=0.85))
dev.off()
