## Create input sceptre object from OTAR input data

# save.image("create_input.rda")
# stop()

# required packages
library(tidyverse)
library(Seurat)
library(sceptre, lib.loc="/home/moonen/R/4.2.2-foss-2022b")

# Process guide metadata ---------------------------------------------------------------------------

# load guide metadata
guide_metadata <- read_csv(snakemake@input$guides, show_col_types = FALSE) %>%
  bind_rows()  # Combine multiple guides files into one dataframe

# get enhancer targeting guide metadata and reformat
enhancer_guides <- guide_metadata %>% 
  filter(grepl(guide_id, pattern = "^chr")) %>% 
  select(grna_id = guide_id, grna_target = peak, chr, start = guide_start, end = guide_end) %>% 
  mutate(grna_id = sub("_", "-", grna_id)) %>% 
  mutate(type = "discovery")

# get non-targeting control guide metadata and reformat
non_targeting_guides <- guide_metadata %>% 
  filter(peak == "Non-targeting") %>% 
  select(grna_id = guide_id, chr, start = guide_start, end = guide_end) %>% 
  mutate(grna_target = "non-targeting", .after = grna_id) %>% 
  mutate(grna_id = sub("_", "-", grna_id)) %>% 
  mutate(type = "non_targeting_control")

# process promoter targeting controls and target id to the name of the target gene
prom_ctrls <- guide_metadata %>%
  filter(grepl(peak, pattern = "_Promoter")) %>% 
  mutate(peak = sub("_Promoter", "", peak)) %>% 
  separate(guide_id, into = c("target", "type", "guide"), sep = "_") %>% 
  mutate(guide_id = paste(type, target, guide, sep = "-")) %>% 
  select(grna_id = guide_id, grna_target = peak, chr, start = guide_start, end = guide_end) %>% 
  mutate(type = "promoter_control")
  
# process enhancer targeting positive controls
enh_ctrls <- guide_metadata %>%
  filter(grepl(guide_id, pattern = "Pilot_positive_control")) %>% 
  mutate(guide_number = sub(".*_([[:digit:]])$", "\\1", guide_id)) %>% 
  mutate(grna_id = paste0("Positive-control-", target_gene, "-Pilot-", guide_number)) %>% 
  select(grna_id, grna_target = target_gene, chr, start = guide_start, end = guide_end) %>% 
  mutate(type = "enhancer_control")

# combine all guide targets information into one table
guide_targets <- bind_rows(enhancer_guides, prom_ctrls, enh_ctrls, non_targeting_guides)

# Extract gene and guide UMI counts ----------------------------------------------------------------

# load merged seurat object containing count data
seurat_obj <- readRDS(snakemake@input$seurat)
seurat_obj <- subset(seurat_obj, doublet_class == "singlet")

# extract UMI counts matrix
counts <- GetAssayData(seurat_obj, assay = "RNA", slot = "counts")

# extract gene and guide counts
guide_pattern <- "^GUIDE-"
guide_rows <- grepl(guide_pattern, rownames(counts))
guide_counts <- counts[guide_rows, ]
gene_counts <- counts[!guide_rows, ]

# remove guide pattern from guide names in guide UMI counts matrix
rownames(guide_counts) <- sub(guide_pattern, "", rownames(guide_counts))

# get guide ids in guide_targets that are also found in guide UMI counts table
guide_ids_in_data <- intersect(rownames(guide_counts), guide_targets$grna_id)
if (length(guide_ids_in_data) != nrow(guide_targets)) {
  warning("Not all guides in metadata found in guide UMI counts table!", call. = FALSE)
}

# extract UMI counts and metadata for guides in both metadata and guide UMI counts table
guide_targets <- filter(guide_targets, grna_id %in% guide_ids_in_data) 
guide_counts <- guide_counts[guide_ids_in_data, ]

# Extract UMI counts for target genes only (optional) ----------------------------------------------

if (snakemake@params$target_genes_only == TRUE) {
  
  # get target genes for this experiment from target genes table
  target_genes <- snakemake@input$target_genes %>% 
    read_tsv(show_col_types = FALSE) %>% 
    filter(panel == snakemake@wildcards$group) %>% 
    pull(gene)
  
  # add all positive control genes (TODO: remove when not needed anymore)
  target_genes <- unique(c(target_genes, prom_ctrls$grna_target, enh_ctrls$grna_target))
  
  # get all target genes found in gene UMI counts
  genes_in_data <- intersect(target_genes, rownames(gene_counts))
  if (all(target_genes %in% genes_in_data) == FALSE) {
    warning("Not all genes in target gene panel found in gene UMI counts table!", call. = FALSE)
  }
  
  # filter gene UMI counts for target genes only
  gene_counts <- gene_counts[genes_in_data, ]
  
}

# Process additional covariates --------------------------------------------------------------------

# load covariates file
covariates <- read_csv(snakemake@input$covars, show_col_types = FALSE,
                       col_types = cols(.default = col_character()))

# Create a cell covariates table for the merged data (adjust if needed)
# cell_covariates <- covariates %>%
#   filter(Name == snakemake@wildcards$group) %>%
#   select(-Notes) %>%
#   mutate(cells = ncol(gene_counts)) %>%
#   uncount(cells) %>%
#   mutate(cell = colnames(gene_counts)) %>%
#   column_to_rownames(var = "cell")

# Assemble sceptre object --------------------------------------------------------------------------

# get all non-mitochondrial genes to use as potential target genes (response)
response_genes <- grep("^MT-", rownames(gene_counts), value = TRUE, invert = TRUE)

# create sceptre object from processed input data
sceptre_object <- import_data(
  response_matrix = gene_counts,
  grna_matrix = guide_counts,
  grna_target_data_frame = guide_targets,
  moi = "low",
  response_names = response_genes
)

# save sceptre object to file
saveRDS(sceptre_object, file = snakemake@output[[1]])