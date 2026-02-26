#####################################################
## Script to select eQTL studies from eQTL catalog ##
#####################################################

#### SETUP ####
#libraries
library(data.table)
library(tidyverse)
library(ggplot2)

#directories
root <- "/g/steinmetz/project/otar_2063/enhancer_screen/"
eqtldir <- paste0(root, "interpretation/eqtl_enrichment/")

# Read data
datasets <- fread(paste0(eqtldir, "data/eQTL_catalogue/tabix_ftp_paths.tsv"))
bulk_studies <- fread(paste0(eqtldir, "data/eQTL_catalogue/Bulk_RNA-seq_studies.tsv"))
sc_studies <- fread(paste0(eqtldir, "data/eQTL_catalogue/Single-cell_RNA-seq_studies.tsv"))

# Select only gene-level results
selected_ftp_links <- datasets %>%
  filter(quant_method == "ge") %>%
  select(ftp_cs_path)

# Select only gene-level results
selected_studies <- datasets %>%
  filter(quant_method == "ge") %>%
  select(study_id, dataset_id, study_label, sample_group, tissue_id, tissue_label, condition_label, sample_size)

#Write list of FTP links
fwrite(selected_ftp_links, paste0(eqtldir, "data/eQTL_catalogue/FTP_links_credible_set_results.tsv"), sep = "\t", col.names = F)
fwrite(selected_studies, paste0(eqtldir, "data/eQTL_catalogue/eQTL_catalogue_studies.tsv"), sep = "\t", col.names = T)
