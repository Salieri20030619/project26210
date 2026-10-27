w <- '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase133_sections2_5_2_7_targeted_support_remote_v1'
d <- read.delim(gzfile(file.path(w,'data/R_CORRELATION_AUDIT_INPUT.tsv.gz')),check.names=FALSE)
a <- read.delim(file.path(w,'data/TYMP_ACTIVITY_PATIENT_EFFECTS.tsv'),check.names=FALSE)
s <- read.delim(file.path(w,'qc/TYMP_ACTIVITY_SUMMARY.tsv'),check.names=FALSE)
maxcor <- 0
for (i in seq_len(nrow(a))) {
 r <- a[i,];g <- d[d$curated_patient_id==r$patient,];x <- rank(g$TYMP_RNA);y <- rank(g[[r$target]])
 if (r$model=='technical_adjusted') {
  c1 <- rank(g$raw_library_size_v1);c2 <- rank(g$raw_detected_genes_v1)
  x <- residuals(lm(x~c1+c2));y <- residuals(lm(y~c1+c2))
 }
 rr <- cor(x,y);maxcor <- max(maxcor,abs(rr-r$effect));stopifnot(abs(rr-r$effect)<1e-9)
}
wp <- function(x) {x<-x[is.finite(x)];x[abs(x)<1e-12]<-0;if(!any(x!=0)) return(1);suppressWarnings(wilcox.test(x,mu=0,exact=FALSE,correct=TRUE)$p.value)}
maxp <- 0
for(i in seq_len(nrow(s))) {
 r<-s[i,];x<-a$effect[a$target==r$target & a$model==r$model];p<-wp(x)
 stopifnot(length(x)==r$n,abs(median(x)-r$median_effect)<1e-10,abs(p-r$P)<1e-10)
 maxp<-max(maxp,abs(p-r$P))
}
for(model in unique(s$model)) {
 ix<-s$type=='metaVIPER_activity' & s$model==model
 stopifnot(max(abs(p.adjust(s$P[ix],'BH')-s$FDR[ix]))<1e-10)
}
cnt<-read.delim(file.path(w,'data/TERMINAL_PSEUDOBULK_COUNTS.tsv'),check.names=FALSE)
b<-read.delim(file.path(w,'data/TERMINAL_ENZYME_PATIENT_EFFECTS.tsv'),check.names=FALSE)
bs<-read.delim(file.path(w,'qc/TERMINAL_ENZYME_SUMMARY.tsv'),check.names=FALSE)
stopifnot(identical(paste(cnt$patient,cnt$method,cnt$gene),paste(b$patient,b$method,b$gene)))
effect<-log2(1+1e6*cnt$late_gene_UMI/cnt$late_total_UMI)-log2(1+1e6*cnt$early_gene_UMI/cnt$early_total_UMI)
stopifnot(max(abs(effect-b$effect))<1e-10)
for(i in seq_len(nrow(bs))) {
 r<-bs[i,];x<-b$effect[b$method==r$method & b$gene==r$gene];p<-wp(x)
 stopifnot(length(x)==r$n,abs(median(x)-r$median_effect)<1e-10,abs(p-r$P)<1e-10)
 maxp<-max(maxp,abs(p-r$P))
}
stopifnot(max(abs(p.adjust(bs$P,'BH')-bs$FDR))<1e-10)
writeLines(c('PASS',paste('Patient correlations checked:',nrow(a)),paste('Terminal count-derived effects checked:',nrow(b)),paste('Max correlation error:',maxcor),paste('Max P error:',maxp),'All primary/sensitivity protein BH families and 39-gene-method terminal BH verified'),file.path(w,'qc/R_INDEPENDENT_VALIDATION.txt'))
capture.output(sessionInfo(),file=file.path(w,'qc/R_SESSION_INFO.txt'))
cat('INDEPENDENT_VALIDATION_PASS\n')
