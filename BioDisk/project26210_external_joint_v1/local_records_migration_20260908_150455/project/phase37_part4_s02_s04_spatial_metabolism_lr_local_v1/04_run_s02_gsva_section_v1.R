suppressPackageStartupMessages({
  library(zellkonverter)
  library(SummarizedExperiment)
  library(Matrix)
  library(GSVA)
  library(BiocParallel)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
getarg <- function(flag) { i <- match(flag, args); if (is.na(i) || i == length(args)) stop("missing ", flag); args[[i+1]] }
manifest_path <- getarg("--manifest"); sample_id <- getarg("--sample"); custom_path <- getarg("--custom-gmt"); kegg_path <- getarg("--kegg-gmt"); out <- getarg("--output")
if (dir.exists(out) || file.exists(out)) stop("refusing to overwrite output: ", out)
dir.create(out, recursive = TRUE)

gmt <- function(path) {
  z <- readLines(path, warn = FALSE, encoding = "UTF-8")
  x <- strsplit(z, "\t", fixed = TRUE)
  ans <- lapply(x, function(v) unique(v[-c(1,2)])); names(ans) <- vapply(x, `[[`, character(1), 1); ans
}
four <- c("Pyrimidine metabolism", "Cysteine and methionine metabolism", "Arginine and proline metabolism", "Glycerophospholipid metabolism")
custom <- gmt(custom_path)[four]; kegg <- gmt(kegg_path)[four]
m <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
r <- m[m$sample == sample_id, , drop = FALSE]; if (nrow(r) != 1L) stop("sample manifest identity failure")

sce <- readH5AD(r$input_h5ad[[1]], use_hdf5 = TRUE, reader = "R")
an <- assayNames(sce); assay_name <- if ("X" %in% an) "X" else an[[1]]
raw <- assay(sce, assay_name)
if (nrow(raw) != r$n_input_genes[[1]] || ncol(raw) != r$n_spots[[1]]) stop("H5AD shape mismatch")
if (is.null(rownames(raw)) || is.null(colnames(raw))) stop("missing gene or spot identifiers")
raw <- as(raw, "dgCMatrix")
if (length(raw@x) && (any(!is.finite(raw@x)) || any(raw@x < 0) || any(abs(raw@x-round(raw@x)) > 1e-8))) stop("X is not raw-count safe")
lib <- Matrix::colSums(raw); if (any(lib <= 0)) stop("zero library size")
expr <- raw
expr@x <- log1p(expr@x * rep(10000/lib, diff(expr@p)))

score_one <- function(sets, source) {
  covered <- lapply(sets, intersect, y = rownames(expr))
  if (any(lengths(covered) < 10L)) stop("pathway coverage below 10")
  param <- ssgseaParam(exprData = expr, geneSets = covered, minSize = 10, maxSize = Inf,
                       alpha = 0.25, normalize = TRUE, checkNA = "auto", use = "everything")
  score <- gsva(param, verbose = FALSE, BPPARAM = SerialParam())
  data.frame(sample = sample_id, spot_id = colnames(score), source = source,
             pathway = rep(rownames(score), each = ncol(score)),
             gsva_ssgsea_score = as.vector(t(score)), stringsAsFactors = FALSE)
}
scores <- rbind(score_one(custom, "CUSTOM_PRIMARY"), score_one(kegg, "KEGG_SENSITIVITY"))
coverage <- do.call(rbind, lapply(list(CUSTOM_PRIMARY=custom,KEGG_SENSITIVITY=kegg), function(sets) do.call(rbind,lapply(names(sets),function(nm) data.frame(sample=sample_id,source=if(identical(sets,custom))"CUSTOM_PRIMARY" else "KEGG_SENSITIVITY",pathway=nm,n_set_genes=length(sets[[nm]]),n_present=sum(sets[[nm]] %in% rownames(expr)))))))

write_atomic <- function(x,path) { part <- paste0(path,".part"); write.table(x,part,sep="\t",quote=FALSE,row.names=FALSE); file.rename(part,path) }
write_atomic(scores,file.path(out,"s02_gsva_spot_scores_v1.tsv")); write_atomic(coverage,file.path(out,"s02_gsva_gene_coverage_v1.tsv"))
payload <- list(status="PASS_S02_GSVA_SECTION",completed_at=format(Sys.time(),"%Y-%m-%dT%H:%M:%S%z"),sample=sample_id,n_spots=ncol(expr),n_genes=nrow(expr),n_score_rows=nrow(scores),method="GSVA 2.4.9 ssgseaParam on raw-count-derived log1p CP10K; alpha=0.25; normalized",figures_created=FALSE)
part <- file.path(out,"status.json.part"); write_json(payload,part,pretty=TRUE,auto_unbox=TRUE); file.rename(part,file.path(out,"status.json")); cat(toJSON(payload,auto_unbox=TRUE,pretty=TRUE),"\n")
