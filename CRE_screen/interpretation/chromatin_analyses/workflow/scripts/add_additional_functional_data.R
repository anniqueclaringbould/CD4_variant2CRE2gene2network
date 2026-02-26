## Add additional functional features to E-G results

# save.image("RDA/add_additional_functional_data.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

## Enhancer-level data -----------------------------------------------------------------------------

# load file containing E-G results with added enhancer chromatin marks
results <- read_csv(snakemake@input$elements, show_col_types = FALSE)

# load eRNA quantifications
erna_counts <- read_tsv(snakemake@input$eRNA_counts, show_col_types = FALSE)

# select eRNA counts for each element
erna_counts_enh <- erna_counts %>% 
  select(grna_target, enh_signal_CPM_NETCAGE = NETCAGE.CPMs,
         enh_signal_TPM_fantomERNA = fantomERNA.TPMs) %>% 
  distinct()

# add eRNA counts to results table
results <- left_join(results, erna_counts_enh, by = "grna_target")

## Gene-level data ---------------------------------------------------------------------------------

# load gene expression, gene constraints and ubiquitous expression tables
gene_tpm <- read_csv(snakemake@input$gene_tpm, show_col_types = FALSE)
gene_constraints <- read_csv(snakemake@input$gene_constraints, show_col_types = FALSE)
ubiq_expression <- read_tsv(snakemake@input$ubiq_expression, show_col_types = FALSE)

# extract relevant data from tpm and ubiquitous expression tables
gene_tpm <- distinct(select(gene_tpm, gene, tpm))
ubiq_expression <- select(ubiq_expression, GeneSymbol, gene_uniq_expr = is_ubiqutious_uniform)

# load gene-level disease information and format columns to add to table
disease_info <- read_tsv(snakemake@input$gene_disease_info, show_col_types = FALSE)
disease_info <- disease_info %>% 
  select(-c(gene_name, chr, gene_start, gene_end)) %>% 
  rename(id = gene_id) %>% 
  setNames(paste0('gene_', names(.)))

# add to results table
results <- results %>% 
  left_join(gene_tpm, by = c("response_id" = "gene")) %>% 
  left_join(gene_constraints, by = "gene_id") %>% 
  left_join(ubiq_expression, by = c("response_id" = "GeneSymbol")) %>% 
  left_join(disease_info, by = "gene_id") %>% 
  rename(gene_tpm = tpm, gene_constraints = lof.pLI)

# save to output file
write_csv(results, snakemake@output[[1]])
