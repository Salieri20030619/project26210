w <- '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase129_sections2_6_2_7_redraw_remote_v1'
q <- file.path(w,'qc')
a <- read.delim(file.path(q,'ALL85_TYMP_PATIENT_DELTAS.tsv'),check.names=FALSE)
s <- read.delim(file.path(q,'ALL85_TYMP_PATIENT_STATISTICS.tsv'),check.names=FALSE)
safe <- function(x) { x <- x[is.finite(x)]; if(length(x)<2) return(NA_real_); if(all(abs(x)<1e-12)) return(1); suppressWarnings(wilcox.test(x,mu=0,exact=FALSE,correct=TRUE)$p.value) }
p <- vapply(s$pathway,function(k) safe(a$delta[a$pathway==k]),numeric(1))
stopifnot(max(abs(p-s$P))<1e-10,max(abs(p.adjust(p,'BH')-s$FDR))<1e-10)
e <- read.delim(file.path(q,'ENZYME_F1_PATIENT_EFFECTS.tsv'),check.names=FALSE)
g <- read.delim(file.path(q,'ENZYME_F1_STATISTICS.tsv'),check.names=FALSE)
pe <- vapply(g$gene,function(k) safe(e$rho[e$gene==k]),numeric(1))
stopifnot(max(abs(pe-g$P))<1e-10,max(abs(p.adjust(pe,'BH')-g$FDR))<1e-10)
writeLines(c('PASS: 85 KEGG paired tests and BH family independently reproduced in R','PASS: 12 enzyme F1 tests and BH family independently reproduced in R',paste('KEGG max P error',max(abs(p-s$P))),paste('Enzyme max P error',max(abs(pe-g$P)))),file.path(q,'R_NEW_STATISTICS_VALIDATION.txt'))
capture.output(sessionInfo(),file=file.path(q,'R_SESSION_INFO.txt'))
cat('R_NEW_STATISTICS_PASS\n')
