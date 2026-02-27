## Create Matrix Market files from count matrices in on-disc seurat objects

# save.image("write_mtx_donor.rda")
# stop()

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(forcats)
  library(Seurat)
  library(BPCells)
  library(Matrix)
  library(rtracklayer)
})

# Load input data ----------------------------------------------------------------------------------

# load seurat object
message("Loading seurat object...")
seu <- readRDS(snakemake@input$seu)

# load genome annotations to extract gene ids for each gene
message("Loading genome annotations...")
annot <- import(snakemake@input$annot, format = "gtf")

# extract gene ids and gene names
genes <- annot[annot$type == "gene"]
gene_ids <- as.data.frame(unique(mcols(genes)[, c("gene_id", "gene_name", "gene_type")]))

# for cases with multiple gene ids per gene name, pick one id
gene_ids <- gene_ids %>% 
  mutate(gene_type = fct_relevel(gene_type, "protein_coding", "lincRNA")) %>% 
  group_by(gene_name) %>% 
  slice_min(order_by = gene_type, with_ties = FALSE) %>% 
  select(-gene_type)

# load list of features to filter out if provided
if (!is.null(snakemake@input$rm_feat)) {
  remove_features <- readLines(snakemake@input$rm_feat)
} else {
  remove_features <- character()
}

# Save counts matrix to file -----------------------------------------------------------------------

# extract data on cells for specified subset (and filtered features if provided)
subset_cells <- grep(paste0("^", snakemake@wildcards$donor, ".+$"), x = colnames(seu))
seu <- subset(seu, cells = subset_cells, features = setdiff(rownames(seu), remove_features))

# convert on-disk counts matrix to sparse matrix
message("Converting counts to in-memory sparse matrix...")
counts <- as(object = seu[["RNA"]]$counts, Class = "dgCMatrix")

# make sure rows in counts matrix are sorted correctly (first genes, then gRNAs)
grna_pattern <- snakemake@params$grna_pattern
gene_rows <- grep(grna_pattern, rownames(counts), invert = TRUE)
grna_rows <- grep(grna_pattern, rownames(counts))
counts <- counts[c(gene_rows, grna_rows), ]

# get temporary mtx output file
tmp_mtx <- tools::file_path_sans_ext(snakemake@output$mtx)

# write count matrix to temp file and compress to create final output file
message("Writing counts to .mtx file...")
invisible(writeMM(counts, file = tmp_mtx))
system2("gzip", args = tmp_mtx)

# Create barcodes and feature files ----------------------------------------------------------------

# create cell barcodes list
message("Creating 'barcodes' and 'features' files...")
barcodes <- data.table(barcode = colnames(counts))

# create features list
features <- data.table(gene_name = rownames(counts))

# add gene ids for genes
features <- left_join(features, gene_ids, by = "gene_name")

# identify all grnas
features <- features %>% 
  mutate(type = if_else(grepl(grna_pattern, gene_name), true = "CRISPR Guide Capture",
                        false = "Gene Expression"))

# set gene id for grnas to grna name and reformat for output
features <- features %>% 
  mutate(gene_id = if_else(type == "CRISPR Guide Capture", true = gene_name, false = gene_id)) %>% 
  select(gene_id, gene_name, type)

# write barcodes and features list to output files
fwrite(barcodes, file = snakemake@output$barcodes, sep = "\t", col.names = FALSE)
fwrite(features, file = snakemake@output$features, sep = "\t", col.names = FALSE)

message("All done!")

# close log file connection
sink()
sink(type = "message")
close(log)
