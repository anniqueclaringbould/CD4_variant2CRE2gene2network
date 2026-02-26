## Extract trans-acting effects from promoter targeting perturbations in enhancer screen

# save.image("RDA/extractTransEnh.rda")
# stop()

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(rtracklayer)
  library(GenomicRanges)
})

# load enhancer screen results
etp <- read_csv(snakemake@input$results, show_col_types = FALSE)

# load genome annotations and extract gene symbols for gene_ids
annot <- import(snakemake@input$annot, format = "gtf")
genes <- as.data.frame(unique(mcols(annot)[, c("gene_id", "gene_name")]))

# extract significant hits involving promoter perturbations
promoter_hits <- etp %>% 
  filter(ctrl_type == "promoter" | element_type == "promoter") %>% 
  filter(significant == TRUE) %>% 
  select(grna_target, overlapping_promoters, response_id) %>% 
  mutate(overlapping_promoters = strsplit(overlapping_promoters, split = ",")) %>% 
  unnest(overlapping_promoters)

# add gene symbol of targeted promoters and only retain trans-acting effects
trans_hits <- promoter_hits %>% 
  left_join(genes, by = c("overlapping_promoters" = "gene_id")) %>% 
  select(promoter = gene_name, trans_gene = response_id) %>% 
  filter(promoter != trans_gene) %>% 
  distinct()
  
# save to output file
write_csv(trans_hits, file = snakemake@output[[1]])
