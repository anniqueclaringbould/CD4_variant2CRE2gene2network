## Create AnnData object with combined RNA counts

# required packages
suppressPackageStartupMessages({
  library(Seurat)
  library(BPCells)
})

# load seruat combined on-disk Seurat object
seu <- readRDS(snakemake@input[[1]])

# create AnnData object file with counts matrix
write_matrix_anndata_hdf5(seu[["RNA"]]$counts, path = snakemake@output[[1]], group = "X",
                          gzip_level = 6)
