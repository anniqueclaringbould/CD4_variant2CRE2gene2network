## Process enhancer screen guide assignments for one panel

# save.image("RDA/process_guide_assignments.rda")
# stop()

suppressPackageStartupMessages({
  library(sceptre)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
})

# load sceptre object containing guide targets data frame
sceptre_object <- readRDS(snakemake@input$sceptre_object)

# extract grna targets table
grna_targets <- sceptre_object@grna_target_data_frame

# load guide assignments matrix
guide_assignments <- readRDS(snakemake@input$guide_assignments)

# set guide assignment matrix column names to cell ids
colnames(guide_assignments) <- rownames(sceptre_object@covariate_data_frame)

# convert guide assignments to long format and only retain information on perturbed cells
guide_assignments <- guide_assignments %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column(var = "grna_id") %>% 
  pivot_longer(cols = -grna_id, names_to = "cell", values_to = "pert") %>% 
  filter(pert == TRUE) %>% 
  select(-pert)

# get targeted element per guide
guide_targets <- sceptre_object@grna_target_data_frame %>% 
  select(grna_id, grna_target) %>% 
  distinct()

# add guide targets to guide assignments
guide_assignments <- guide_assignments %>% 
  left_join(guide_targets, by = "grna_id") %>% 
  select(cell, grna_id, grna_target)

# add number of guides per cell
guide_assignments <- add_count(guide_assignments, cell, name = "guides_per_cell")

# compute guide assignment summary statistics
summary_stats <- tibble(
  total_cells = n_distinct(rownames(sceptre_object@covariate_data_frame)),
  cells_with_guide = n_distinct(guide_assignments$cell),
  cells_with_1_guide = nrow(filter(guide_assignments, guides_per_cell == 1)),
  pct_cells_with_guide = cells_with_guide / total_cells,
  pct_cells_with_1_guide = cells_with_1_guide / total_cells
)

# write all output to giles files
write_tsv(guide_assignments, file = snakemake@output$guide_assignments)
write_tsv(summary_stats, file = snakemake@output$summary_stats)
write_tsv(grna_targets, file = snakemake@output$grna_targets)
