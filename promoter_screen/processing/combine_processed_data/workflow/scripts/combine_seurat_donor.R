## Merge seurat objects for one donor

# save.image("combine_seurat.rda")
# stop()

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(BPCells)
})

# function to load seurat object and create minimal new object containing only the raw counts
load_minimal_seurat <- function(x, name) {
  message("\tLoading ", x)
  seu <- readRDS(x)
  seu_minimal <- CreateSeuratObject(counts = GetAssayData(seu, assay = "RNA", layer = "counts"),
                                    project = name)
  return(seu_minimal)
}

# get all seurat objects to merge
seurat_objects <- unlist(snakemake@input)
names(seurat_objects) <- basename(dirname(seurat_objects))

# load all seurat objects
message("Loading input files...")
input_objects <- lapply(names(seurat_objects), FUN = function(x) {
  load_minimal_seurat(seurat_objects[[x]], name = x)
})
names(input_objects) <- names(seurat_objects)

# free up unused memory
invisible(gc())

# merge seurat objects into one
message("Merging seurat objects...")
merged_object <- merge(
  x = input_objects[[1]],
  y = input_objects[-1],
  add.cell.ids = names(input_objects),
  project = "promoter_screen"
)

# combine RNA layer across donors
merged_object[["RNA"]] <- JoinLayers(merged_object[["RNA"]])

# make sure the data type of the counts matrix is integer 
merged_object[["RNA"]]$counts <- convert_matrix_type(merged_object[['RNA']]$counts, "uint32_t")

# write counts matrix to on-disk directory and replace matrix in seurat object with on-disk matrix
write_matrix_dir(mat = merged_object[["RNA"]]$counts, dir = snakemake@output$counts_dir)
counts_mat <- open_matrix_dir(dir = snakemake@output$counts_dir)
merged_object[["RNA"]]$counts <- counts_mat

# save merged seurat object to file
message("Writing to output files...")
saveRDS(merged_object, file = snakemake@output[[1]])

message("All done!")

# close log file connection
sink()
sink(type = "message")
close(log)
