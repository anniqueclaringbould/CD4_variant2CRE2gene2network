## Combine Seurat objects from donors

# save.image("combine_all_donors.rda")
# stop()

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
})

# load all seurat objects using on-disk count matrices to merge
message("Loading input files...")
seurat_objects <- unlist(snakemake@input)
input_objects <- lapply(seurat_objects, FUN = readRDS)

# convert all on-disk matrices to "double" for merging
input_objects <- lapply(input_objects, FUN = function(x) {
  x[["RNA"]]$counts <- convert_matrix_type(x[["RNA"]]$counts, "double")
  return(x)
})

# merge all seurat objects into one
message("Merging seurat objects...")
merged_object <- merge(x = input_objects[[1]], y = input_objects[-1], project = "promoter_screen")

message("Combining RNA layer...")
merged_object[["RNA"]] <- JoinLayers(merged_object[["RNA"]])

# write combined matrix to new on-disk files and link to merged seurat object
message("Writing merged counts matrix to on-disk files...")
merged_object[["RNA"]]$counts <- convert_matrix_type(merged_object[['RNA']]$counts, "uint32_t")
write_matrix_dir(mat = merged_object[["RNA"]]$counts, dir = snakemake@output$counts_dir)
counts_mat <- open_matrix_dir(dir = snakemake@output$counts_dir)
merged_object[["RNA"]]$counts <- counts_mat

# save merged Seurat object to file
message("Writing combined Seurat object to output file...")
saveRDS(merged_object, file = snakemake@output$seu)

message("All done!")

# close log file connection
sink()
sink(type = "message")
close(log)
