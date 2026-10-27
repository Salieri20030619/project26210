options(stringsAsFactors=FALSE)
P <- '/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1/phase295_literature_guided_sequential_trials_20260930_v1'
O <- file.path(P,'step1_receiver_network')
gate <- function(){stopifnot(system(paste('/mnt/BioDisk/ludingyu/project26210/check_memory_gate.sh >>',shQuote(file.path(P,'logs/memory_gate.log')),'2>&1'))==0)}
gate()
suppressPackageStartupMessages(library(nichenetr))
r <- '/mnt/BioDisk/ludingyu/project26210/external/nichenet'
lt <- readRDS(file.path(r,'ligand_target_matrix_nsga2r_final.rds'))
lr <- unique(as.data.frame(readRDS(file.path(r,'lr_network_human_21122021.rds')))[,c('from','to')])
d <- read.delim(file.path(O,'RECEIVER_GENE_INPUT.tsv'))
s <- read.delim(file.path(O,'SENDER_GENE_INPUT.tsv'))
flag <- function(x)tolower(as.character(x))=='true'
bg <- intersect(d$gene[flag(d$background)],rownames(lt))
gs <- intersect(d$gene[flag(d$target)],bg)
rec <- intersect(d$gene[flag(d$receiver_expressed)],lr$to)
lig <- intersect(unique(lr$from[lr$to %in% rec]),colnames(lt))
focus <- intersect(lig,s$gene[flag(s$macrophage_expressed)])
stopifnot(length(bg)>100,length(gs)>=10,length(gs)<length(bg),length(lig)>5)
writeLines(bg,file.path(O,'BACKGROUND_FROZEN.txt'));writeLines(gs,file.path(O,'TARGETS_FROZEN.txt'))
writeLines(lig,file.path(O,'ALL_POTENTIAL_LIGANDS_FROZEN.txt'))
writeLines(focus,file.path(O,'MACROPHAGE_LIGANDS_FROZEN.txt'))
write.table(data.frame(targets=nrow(lt),ligands=ncol(lt),background=length(bg),geneset=length(gs),eligible_ligands=length(lig),focused_ligands=length(focus),nichenetr=as.character(packageVersion('nichenetr'))),file.path(O,'NETWORK_DIMENSIONS.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
gate()
a <- as.data.frame(predict_ligand_activities(geneset=gs,background_expressed_genes=bg,ligand_target_matrix=lt,potential_ligands=lig))
a <- a[order(-a$pearson,a$test_ligand),]
a$rank_all <- seq_len(nrow(a));a$macrophage_expressed <- a$test_ligand %in% focus
a$c4_expressed <- a$test_ligand %in% s$gene[flag(s$c4_expressed)]
a$rank_macrophage <- NA_integer_;a$rank_macrophage[a$macrophage_expressed] <- seq_len(sum(a$macrophage_expressed))
write.table(a,file.path(O,'LIGAND_ACTIVITIES_ALL.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
pairs <- lr[lr$from %in% lig & lr$to %in% rec,];write.table(pairs,file.path(O,'EXPRESSED_RECEPTOR_LINKS.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
chosen <- unique(c(head(a$test_ligand[a$macrophage_expressed],20),intersect(c('SPP1','TGFB1','IL1A','IL1B'),lig)))
links <- do.call(rbind,lapply(chosen,function(l){v <- lt[gs,l];j <- order(-v,names(v));data.frame(ligand=l,target=names(v)[j],prior_regulatory_potential=as.numeric(v[j]),target_rank=seq_along(j))}))
write.table(links,file.path(O,'CANDIDATE_TARGET_LINKS.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
# Independent Pearson recomputation: binary actual target membership against each prior column.
truth <- as.numeric(bg %in% gs)
manual <- vapply(a$test_ligand,function(l)cor(truth,lt[bg,l]),numeric(1))
stopifnot(max(abs(manual-a$pearson),na.rm=TRUE)<1e-10)
stopifnot(!any(gs %in% d$gene[flag(d$excluded_definition_gene)]),all(gs %in% bg),!anyDuplicated(a$test_ligand),setequal(a$test_ligand,lig))
write.table(data.frame(check=c('manual_Pearson_max_error','no_definition_genes','target_subset_background','all_eligible_ligands_ranked'),value=c(max(abs(manual-a$pearson),na.rm=TRUE),TRUE,TRUE,TRUE)),file.path(O,'NETWORK_VALIDATION.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(O,'R_SESSION_INFO.txt'))
cat('STEP1_NETWORK_COMPLETE\n');print(head(a,10));print(a[a$test_ligand %in% c('SPP1','TGFB1','IL1A','IL1B'),])
