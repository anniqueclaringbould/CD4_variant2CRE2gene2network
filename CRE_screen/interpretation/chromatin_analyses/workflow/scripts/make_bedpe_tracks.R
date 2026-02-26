## Create .bedpe tracks from CRISPR results for visualization in IGV genome browser

# save.image("RDA/tracks.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

# Load and filter ETPs -----------------------------------------------------------------------------

# load ETP results table
results <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

# for control perturbations which were in multiple panels, select the most significant one
controls <- results %>% 
  filter(ctrl_perturbation == TRUE) %>% 
  group_by(response_id, grna_target) %>% 
  slice_min(p_value, n = 1, with_ties = FALSE)

# add back to results table with discovery pairs
results <- results %>% 
  filter(ctrl_perturbation == FALSE) %>% 
  bind_rows(controls) %>% 
  arrange(desc(p_value))

# filter out any trans-chromosomal interactions
results <- filter(results, pert_chr == gene_chr)

# Make elements track ------------------------------------------------------------------------------

# extract unique elements and label each one based on if it occurred in any significant
# enhancer-like interactions
elements <- results %>% 
  mutate(color = "195,100,92") %>% 
  select(pert_chr, pert_start, pert_end, grna_target, enhancer_like_interaction, color) %>% 
  distinct() %>% 
  group_by(grna_target) %>% 
  slice_max(order_by = enhancer_like_interaction, with_ties = FALSE)

# create a track with all crispr elements
elements_track <- elements %>% 
  mutate(score = 0, strand = ".", .after = grna_target) %>%
  select(chrom = pert_chr, chromStart = pert_start, chromEnd = pert_end, name = grna_target, 
         score, strand, thickStart = pert_start, thickEnd = pert_end, itemRgb = color)

# save to bed file
elements_header <- 'track name="CrisprElements" itemRgb="On"'
writeLines(elements_header, con = snakemake@output$elements)
write_tsv(elements_track, snakemake@output$elements, append = TRUE)

# Make E-G pairs track -----------------------------------------------------------------------------

# define color based on whether a given pair is significant and effect size direction 
results <- results %>% 
  mutate(color = case_when(
    enhancer_like_interaction == FALSE ~ "190,190,190",
    log_2_fold_change <= 0 ~ "7,186,185",
    log_2_fold_change > 0 ~ "231,139,36",
    TRUE ~ "0,0,0"
  ))

# create bedpe track table for significant pairs
bedpe <- results %>% 
  mutate(name = paste0(response_id, "|", grna_target), pert_strand = ".") %>% 
  select(chr1 = pert_chr, x1 = pert_start, x2 = pert_end, chrom2 = gene_chr, start2 = gene_tss,
         end2 = gene_tss, name, score = log_2_fold_change, strand1 = pert_strand,
         strand2 = gene_strand, color) %>% 
  mutate(start2 = start2 - 1)

# save tracks to output files
interactions_header <- "chr1	x1	x2	chrom2	start2	end2	name	score	strand1	strand2	color"
writeLines(interactions_header, con = snakemake@output$sig_bedpe)
write_tsv(filter(bedpe, color != "190,190,190"), file = snakemake@output$sig_bedpe, append = TRUE)

writeLines(interactions_header, con = snakemake@output$nonsig_bedpe)
write_tsv(filter(bedpe, color == "190,190,190"), file = snakemake@output$nonsig_bedpe,
          append = TRUE)
