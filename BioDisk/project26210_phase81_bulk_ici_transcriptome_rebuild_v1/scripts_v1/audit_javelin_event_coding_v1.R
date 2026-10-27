#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
})

args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 2L) stop("Usage: audit_javelin_event_coding_v1.R clinical.tsv output.tsv")
clinical <- fread(args[[1L]])
rows <- list()
index <- 0L
for (arm in c("Avelumab+Axitinib", "Sunitinib")) {
  data <- clinical[TRT01P == arm & is.finite(PFS_P) & PFS_P > 0 & PFS_P_CNSR %in% c(0, 1)]
  for (coding in c("DIRECT_CNSR_AS_EVENT", "INVERTED_ONE_MINUS_CNSR")) {
    event <- if (coding == "DIRECT_CNSR_AS_EVENT") data$PFS_P_CNSR else 1 - data$PFS_P_CNSR
    fit <- survfit(Surv(data$PFS_P, event) ~ 1)
    index <- index + 1L
    rows[[index]] <- data.table(
      arm=arm,
      coding=coding,
      n=nrow(data),
      n_event=sum(event),
      median_pfs_months=unname(summary(fit)$table[["median"]]),
      observed_status_ones=sum(data$PFS_P_CNSR == 1),
      observed_status_zeros=sum(data$PFS_P_CNSR == 0)
    )
  }
}
result <- rbindlist(rows)
fwrite(result, args[[2L]], sep="\t")
print(result)
