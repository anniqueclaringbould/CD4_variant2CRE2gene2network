## Create table with additional cell covariates for sceptre

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Seurat)
})

# load seurat object
seu <- readRDS(snakemake@input$seu)

# extract cell barcodes and all additionally specified covariates from seurat object meta data
covariates <- seu@meta.data %>% 
  rownames_to_column(var = "cell_barcode") %>% 
  select(all_of(c("cell_barcode", unlist(snakemake@params$seurat_covariates))))

# infer batch cell covariates from cell barcodes
covariates <- covariates %>% 
  separate(cell_barcode, into = c("donor_chip", "panel", "cell"), remove = FALSE, sep = "_") %>% 
  separate(donor_chip, into = c("donor", "chip"), sep = 1) %>% 
  select(-cell)

# if specified add, PCA cell embeddings to covariate data frame
if (!is.null(snakemake@params$seurat_pcs)) {
  pc_covars <- seu@reductions$pca.full@cell.embeddings[, snakemake@params$seurat_pcs] %>% 
    as.data.frame() %>% 
    rownames_to_column(var = "cell_barcode")
  covariates <- left_join(covariates, pc_covars, by = "cell_barcode")
}

# add additional covariates if specified and add to batch covariates table
if (!is.null(snakemake@input$add_covars)) {
  add_covars <- readRDS(snakemake@input$add_covars)
  covariates <- left_join(covariates, add_covars, by = "cell_barcode")
}

# convert cell barcodes to row names for sceptre
covariates <- column_to_rownames(covariates, var = "cell_barcode")
  
# save cell covariates table to file
saveRDS(covariates, file = snakemake@output[[1]])
