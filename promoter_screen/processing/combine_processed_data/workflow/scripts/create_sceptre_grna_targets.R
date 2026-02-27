## Create gRNA targets table for sceptre

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# load all feature files and combine into one table (these are usually all identical)
feature_cols <- c("gene_id", "gene_symbol", "type")
features <- snakemake@input %>% 
  lapply(FUN = read_tsv, col_names = feature_cols, show_col_types = FALSE) %>% 
  bind_rows() %>% 
  distinct()

# extract all gRNAs based on the type column in the features list
grnas <- filter(features, type == "CRISPR Guide Capture")

# infer gRNA target genes from guide id
grna_targets <- grnas %>% 
  select(grna_id = gene_id) %>% 
  mutate(grna_target = sub("GUIDE-(.+)-[[:digit:]]+", "\\1", grna_id))

# set correct gRNA target for non-targeting control guides
grna_targets <- grna_targets %>% 
  mutate(grna_target = if_else(grna_target == "NO-TARGET", true = "non-targeting",
                               false = grna_target))

# save gRNA targets table to file
write_tsv(grna_targets, file = snakemake@output[[1]])
