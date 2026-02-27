## Make tailored locus network plot for TYK2-CDC37 and DEXI loci

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(here)
})

# locus network functions
source("locus_network_functions.R")

# Required input files -----------------------------------------------------------------------------

# file containing the annotated CRE screen results. this is produced by the 
# 'CRE_screen/interpretation/chromatin_analyses' workflow
enh_results_file <- "results_df_with_promoterC_annotated_allFeatures.csv"

# file containing all promoter screen hits. this is produced by running the
# `1_get_promoter_screen_hits.R` to extract significant hits from DESeq2 outputs produced by code in
# `promoter_screen/interpretation/DEG_testing`
prom_hits_file <- here("promoter_screen_hits.tsv.gz")

# Process input data -------------------------------------------------------------------------------

# load enhancer screen results
enh_results <- fread(enh_results_file)

# filter for significant enhancer-like interactions and reformat for network functions
chr_levels <- paste0("chr", c(as.character(seq(1:22)), "MT", "X", "Y"))
enh_hits <- enh_results %>%
  filter(enhancer_like_interaction == TRUE) %>% 
  select(regulator = grna_target, target = response_id, effect_size = log_2_fold_change,
         pert_chr, pert_start, pert_end, gene_chr, gene_tss, dist_to_tss) %>% 
  mutate(pert_chr = factor(pert_chr, levels = chr_levels, ordered = TRUE),
         gene_chr = factor(pert_chr, levels = chr_levels, ordered = TRUE))

# fix an issue where the same gene has multiple TSS annotations...
gene_tss <- enh_hits %>% 
  select(target, gene_tss) %>% 
  distinct() %>% 
  group_by(target) %>% 
  slice_head(n = 1)

# add unique TSS coordinates to enh_hits
enh_hits <- enh_hits %>% 
  select(-gene_tss) %>% 
  left_join(gene_tss, by = "target") %>% 
  relocate(gene_tss, .after = "gene_chr")

# load promoter screen results
prom_hits <- fread(prom_hits_file)

# only retain hits with at least an LFC of 0.2 and reformat for network functions
prom_hits <- prom_hits %>% 
  filter(!contrast == variable) %>% 
  filter(abs(log_fc) >= 0.2) %>% 
  select(regulator = contrast, target = variable, effect_size = log_fc) %>% 
  mutate(pert_chr = NA_character_, pert_start = NA_integer_, pert_end = NA_integer_,
         gene_chr = NA_character_, gene_tss = NA_integer_)

## Make locus plot ---------------------------------------------------------------------------------

# Make locus plots and save them to pdfs. Due to how pdfs are written the legend is hidden by
# clipping masks, which can be removed by editing the files in e.g. Adobe Illustrator

# TYK2-CD37 locus
pdf("tyk2_cd37_locus_plot.pdf", width = 10, height = 8)
make_locus_network(enh_hits, prom_hits, start_gene = c("TYK2", "CDC37"), layers = 2,
                   node_size = 3, cluster_gap = 0.4, min_node_spacing = 1, legend_cex = 0.6,
                   jitter_amount = 0)
dev.off()


# DEXI locus
pdf("dexi_locus_plot.pdf", width = 10, height = 8)
make_locus_network(enh_hits, prom_hits, start_gene = "DEXI", layers = 2, node_size = 3,
                   cluster_gap = 0.4, min_node_spacing = 1, legend_cex = 0.6, jitter_amount = 0)
dev.off()
