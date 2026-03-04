library(Seurat)
library(dplyr)
library(scDblFinder)
library(SingleCellExperiment)  # Required for scDblFinder
library(data.table)  # For fwrite

# Snakemake input variables
seurat_files <- snakemake@input$seurat_files
covars <- snakemake@input$covars
target_genes <- snakemake@input$target_genes
merged_output_file <- snakemake@output$merged_seurat
doublet_output_file <- snakemake@output$doublets_filtered
group <- snakemake@wildcards$group
doublet_rate <- snakemake@params$doublet_rate  # Corrected access to doublet rate

# Read the covariates and target genes files
covariates <- read.csv(covars)
target_gene_panels <- read.table(target_genes, header = TRUE, sep = "\t")

# Print the column names to verify structure
print(colnames(target_gene_panels))

# Check if 'panel' column exists in target_gene_panels
if(!"panel" %in% colnames(target_gene_panels)) {
  stop("The 'panel' column is not found in target_gene_panels.")
}

# Filter target genes for the current group
target_genes_for_group <- target_gene_panels %>%
  filter(panel == group) %>%
  pull(gene)

# Initialize an empty list to store Seurat objects
datasets <- list()

# Loop through the Seurat files and read them
for (i in seq_along(seurat_files)) {
  file_path <- seurat_files[i]
  seurat_obj <- readRDS(file_path)
  sample_name <- gsub("^Enhancerscreen-2-2-|_Seurat.rds$", "", basename(file_path))
  seurat_obj@meta.data$orig.ident <- sample_name
  
# Calculate doublet rate and find doublets
  multiplet_rate <- doublet_rate[i]
  # Find doublets with this rate
  sce <- scDblFinder(SingleCellExperiment(list(counts=GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))),
                     samples=seurat_obj@meta.data$orig.ident,
                     dbr=multiplet_rate/100)
  
  seurat_obj$doublet_class <- sce$scDblFinder.class
  seurat_obj$doublet_score <- sce$scDblFinder.score
  
  datasets[[sample_name]] <- seurat_obj
}

# Merge the Seurat objects
seur <- datasets[[1]]
for (i in 2:length(datasets)) {
  seur <- merge(seur, y = datasets[[i]], add.cell.ids = c(names(datasets)[1], names(datasets)[i]))
}

# Save the merged Seurat object
saveRDS(seur, file = merged_output_file)

# Save the doublet information
doublet_info <- as.data.frame(table(seur$doublet_class))
fwrite(doublet_info, file = doublet_output_file)