#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
  library(dorothea)
})
PROJECT <- "/mnt/BioDisk/wangmingyuan/project26210_external_joint_v1_phase64_large_immune_validation_v1"
OUT <- file.path(PROJECT, "inputs_v1/GSE207493/tf_regulon_de_novo_v1/public_regulons_v1")
dir.create(OUT, recursive=TRUE, showWarnings=FALSE)
proxy <- "http://127.0.0.1:7897"
Sys.setenv(http_proxy=proxy, https_proxy=proxy, HTTP_PROXY=proxy, HTTPS_PROXY=proxy, ALL_PROXY=proxy)

collectri <- as.data.table(decoupleR::get_collectri(organism="human", split_complexes=FALSE))
if (any(!c("source", "target", "mor") %in% names(collectri))) stop("Unexpected CollecTRI schema")
if (!"likelihood" %in% names(collectri)) collectri[, likelihood := 1.0]
collectri <- unique(collectri[, .(source=as.character(source), target=as.character(target),
                                  mor=as.numeric(mor), likelihood=as.numeric(likelihood))])
if (nrow(collectri) < 30000L || uniqueN(collectri$source) < 500L) stop("CollecTRI unexpectedly small")
fwrite(collectri, file.path(OUT, "collectri_human_frozen_v1.tsv.gz"), sep="\t")

data(dorothea_hs)
dorothea_abc <- as.data.table(dorothea_hs)[confidence %chin% c("A", "B", "C")]
dorothea_abc[, likelihood := fifelse(confidence == "A", 1.0,
                              fifelse(confidence == "B", 0.5, 1/3))]
dorothea_abc <- unique(dorothea_abc[, .(source=as.character(tf), target=as.character(target),
  mor=as.numeric(mor), likelihood=as.numeric(likelihood), confidence=as.character(confidence))])
if (nrow(dorothea_abc) < 10000L || uniqueN(dorothea_abc$source) < 100L) stop("DoRothEA ABC unexpectedly small")
fwrite(dorothea_abc, file.path(OUT, "dorothea_human_ABC_frozen_v1.tsv.gz"), sep="\t")

summary <- data.table(
  resource=c("CollecTRI", "DoRothEA_ABC"), n_edges=c(nrow(collectri), nrow(dorothea_abc)),
  n_sources=c(uniqueN(collectri$source), uniqueN(dorothea_abc$source)),
  n_targets=c(uniqueN(collectri$target), uniqueN(dorothea_abc$target)),
  source_method=c("decoupleR::get_collectri human split_complexes=FALSE via proxy",
                  "dorothea::dorothea_hs confidence A/B/C packaged data")
)
fwrite(summary, file.path(OUT, "public_regulon_summary_v1.tsv"), sep="\t")
writeLines(capture.output(sessionInfo()), file.path(OUT, "session_info_v1.txt"))
writeLines(format(Sys.time(), tz="", usetz=TRUE), file.path(OUT, "RUN_COMPLETE"))
print(summary)
cat("PUBLIC_REGULONS_V1_PASS\n")
