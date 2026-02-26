## Annotate enhancer screen results with likely trans hits based on perturbation targeting annotated
## promoters, have effects on genes on other chromosomes or overlap genes, for which promoter
## perturbations have a significant effect on the same target genes

# save.image("RDA/annotate_trans_hits_enhancer_screen.rda")
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

# load trans-acting hits from promoter targeting perturbations
trans_hits_enh  <- read_csv(snakemake@input$enh_hits,  show_col_types = FALSE)
trans_hits_prom <- read_csv(snakemake@input$prom_hits, show_col_types = FALSE)
trans_hits <- distinct(bind_rows(enhancer = trans_hits_enh, promoter = trans_hits_prom))

# load genome annotations and extract gene symbols for gene_ids
annot <- import(snakemake@input$annot, format = "gtf")
genes <- as.data.frame(unique(mcols(annot)[, c("gene_id", "gene_name")]))

# get all promoters leading to trans effects for each gene
trans_hits <- trans_hits %>% 
  group_by(trans_gene) %>% 
  summarize(trans_effect_promoters = list(promoter), .groups = "drop")

# add gene symbols of genes overlapping perturbed elements
overlapping_genes <- etp %>% 
  select(response_id, grna_target, panel, overlapping_genes) %>% 
  mutate(overlapping_genes = strsplit(overlapping_genes, split = ",")) %>% 
  unnest(overlapping_genes) %>% 
  left_join(genes, by = c("overlapping_genes" = "gene_id"))
  
# collapse into list
overlapping_genes <- overlapping_genes %>% 
  group_by(response_id, grna_target, panel) %>% 
  summarize(overlapping_gene_symbols = list(gene_name), .groups = "drop")

# add back to results table
etp <- etp %>% 
  left_join(overlapping_genes, by = c("response_id", "grna_target", "panel")) %>% 
  relocate(overlapping_gene_symbols, .after = overlapping_genes)

# add promoters of genes involved in trans hits of target gene of each ETP
etp <- left_join(etp, trans_hits, by = c("response_id" = "trans_gene"))

# annotate likely trans effects based on element type and location, or if perturbed element overlaps
# a gene that showed a trans effect on the same target gene
etp <- etp %>%
  rowwise() %>% 
  mutate(likely_trans_effect = case_when(
    pert_chr != gene_chr ~ TRUE,
    element_type == "promoter" ~ TRUE,
    dist_to_gene == 0 ~ TRUE,
    any(trans_effect_promoters %in% overlapping_gene_symbols) ~ TRUE,
    TRUE ~ FALSE
  ))

# save annotated table to output file
write_csv(etp, file = snakemake@output[[1]])
