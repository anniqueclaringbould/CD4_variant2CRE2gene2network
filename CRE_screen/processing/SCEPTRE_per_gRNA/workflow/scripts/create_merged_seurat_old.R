library(Seurat)
library(dplyr)

# Snakemake input variables
seurat_files <- snakemake@input$seurat_files
covars <- snakemake@input$covars
target_genes <- snakemake@input$target_genes
output_file <- snakemake@output[[1]]
group <- snakemake@wildcards$group  # Ensure this is correctly passed

# Read the covariates and target genes files
covariates <- read.csv(covars)
target_gene_panels <- read.table(target_genes, header = TRUE, sep = "\t")

# Print the column names to verify structure
print(colnames(target_gene_panels))

# Check if 'panel' column exists in target_gene_panels
if(!"panel" %in% colnames(target_gene_panels)) {
  stop("The 'panel' column is not found in target_gene_panels.")
}

# Filter target genes for the current group (using 'panel' instead of 'group')
target_genes_for_group <- target_gene_panels %>%
  filter(panel == group) %>%
  pull(gene)

# Initialize an empty list to store Seurat objects
datasets <- list()

# Loop through the Seurat files and read them
for (file_path in seurat_files) {
  seurat_obj <- readRDS(file_path)
  sample_name <- gsub("^Enhancerscreen-2-2-|_Seurat.rds$", "", basename(file_path))
  seurat_obj@meta.data$orig.ident <- sample_name
  datasets[[sample_name]] <- seurat_obj
}

# Merge the Seurat objects
seur <- datasets[[1]]
for (i in 2:length(datasets)) {
  seur <- merge(seur, y = datasets[[i]], add.cell.ids = c(names(datasets)[1], names(datasets)[i]))
}

# Save the merged Seurat object
saveRDS(seur, file = output_file)