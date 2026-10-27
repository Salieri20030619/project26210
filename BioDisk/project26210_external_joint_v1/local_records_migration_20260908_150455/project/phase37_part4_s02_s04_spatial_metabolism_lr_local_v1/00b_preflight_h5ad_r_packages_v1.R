packages <- c("zellkonverter", "HDF5Array", "rhdf5", "reticulate", "Matrix")
out <- data.frame(package = packages, available = vapply(packages, requireNamespace, logical(1), quietly = TRUE))
write.table(out, file = stdout(), sep = "\t", quote = FALSE, row.names = FALSE)
for (p in packages[out$available]) cat(p, as.character(packageVersion(p)), "\n", sep = "\t")
