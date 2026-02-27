## Create positive control pairs based on gRNAs targeting a gene promoter

# save.image("positive_controls.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

# load gRNA targets
grna_targets <- read_tsv(snakemake@input$grna_targets, show_col_types = FALSE)

# load all feature files and combine into one table (these are usually all identical)
feature_cols <- c("gene_id", "gene_symbol", "type")
features <- snakemake@input$features %>% 
  lapply(FUN = read_tsv, col_names = feature_cols, show_col_types = FALSE) %>% 
  bind_rows() %>% 
  distinct()

# get gene_ids for all genes
gene_ids <- features %>% 
  filter(type == "Gene Expression") %>% 
  select(gene_symbol, gene_id)

# add gene ids to gRNA targets
grna_targets <- left_join(grna_targets, gene_ids, by = c("grna_target" = "gene_symbol"))

# create positive controls by pairing promoter perturbations with they're own target gene
pos_ctrls <- grna_targets %>% 
  select(grna_target, response_id = gene_id) %>% 
  distinct() %>% 
  drop_na()

# save positive control pairs to file
write_tsv(pos_ctrls, file = snakemake@output[[1]])
