#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

features_1 <- c("a", "b", "c")
matrix_1 <- sparseMatrix(i = c(1, 3), j = c(1, 2), x = c(1, 2), dims = c(3, 2),
                         dimnames = list(features_1, c("x1", "x2")))
features_2 <- c("c", "a")
matrix_2 <- sparseMatrix(i = c(1, 2), j = c(1, 2), x = c(3, 4), dims = c(2, 2),
                         dimnames = list(features_2, c("y1", "y2")))
chunks <- lapply(list(matrix_1, matrix_2), function(m) {
  list(features = rownames(m), i = m@i, p = m@p, x = m@x, nrow = nrow(m), ncol = ncol(m))
})
feature_union <- unique(unlist(lapply(chunks, `[[`, "features"), use.names = FALSE))
total_nnz <- sum(vapply(chunks, function(z) length(z$x), integer(1)))
total_cells <- sum(vapply(chunks, function(z) z$ncol, integer(1)))
x_out <- numeric(total_nnz)
i_out <- integer(total_nnz)
p_out <- integer(total_cells + 1)
nnz_offset <- 0L
col_offset <- 0L

for (chunk in chunks) {
  row_map <- match(chunk$features, feature_union)
  chunk_i <- chunk$i
  chunk_p <- chunk$p
  chunk_x <- chunk$x
  remapped_needed <- FALSE
  if (length(row_map) > 1L && any(diff(row_map) <= 0L)) {
    old_matrix <- methods::new("dgCMatrix", i = chunk$i, p = chunk$p, x = chunk$x,
      Dim = as.integer(c(chunk$nrow, chunk$ncol)), Dimnames = list(chunk$features, NULL))
    row_projection <- sparseMatrix(i = row_map, j = seq_along(row_map), x = 1,
      dims = c(length(feature_union), length(row_map)))
    remapped <- as(row_projection %*% old_matrix, "dgCMatrix")
    chunk_i <- remapped@i
    chunk_p <- remapped@p
    chunk_x <- remapped@x
    remapped_needed <- TRUE
  }
  n_chunk_nnz <- length(chunk_x)
  idx <- seq.int(nnz_offset + 1, nnz_offset + n_chunk_nnz)
  x_out[idx] <- chunk_x
  if (!remapped_needed) i_out[idx] <- row_map[chunk_i + 1L] - 1L else i_out[idx] <- chunk_i
  p_idx <- seq.int(col_offset + 2L, col_offset + chunk$ncol + 1L)
  p_out[p_idx] <- as.integer(chunk_p[-1L] + nnz_offset)
  nnz_offset <- as.integer(nnz_offset + n_chunk_nnz)
  col_offset <- col_offset + chunk$ncol
}

assembled <- methods::new("dgCMatrix", i = i_out, p = p_out, x = x_out,
  Dim = as.integer(c(length(feature_union), total_cells)),
  Dimnames = list(feature_union, c(colnames(matrix_1), colnames(matrix_2))))
validObject(assembled)
expected_2 <- sparseMatrix(
  i = match(rownames(matrix_2), feature_union)[matrix_2@i + 1L],
  j = rep.int(seq_len(ncol(matrix_2)), diff(matrix_2@p)), x = matrix_2@x,
  dims = c(length(feature_union), ncol(matrix_2)), dimnames = list(feature_union, colnames(matrix_2))
)
expected <- cbind(matrix_1, expected_2)
stopifnot(identical(as.matrix(assembled), as.matrix(expected)), is.integer(assembled@p))
cat("SPARSE_ASSEMBLY_UNIT_TEST=PASS\n")
cat("P_SLOT_TYPE=", typeof(assembled@p), "\n", sep = "")
cat("FEATURE_UNION=", paste(feature_union, collapse = ","), "\n", sep = "")
