## Calculate miscellaneous enhancer screen numbers

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# load guide assignments and summary statistics files
guide_assignments <- read_tsv(snakemake@input$guide_assignments, show_col_types = FALSE)
summary_stats <- read_tsv(snakemake@input$summary_stats, show_col_types = FALSE)

# calculate average of summary stats across all panels
summary_stats_across_panels <- summary_stats %>% 
  summarize(
    total_cells = sum(total_cells),
    cells_with_guide = sum(cells_with_guide),
    cells_with_1_guide = sum(cells_with_1_guide),
    pct_cells_with_guide = cells_with_guide / total_cells,
    pct_cells_with_1_guide = cells_with_1_guide / total_cells
    )

# filter guide assignments to cells with exactly one assigned guide
guide_assignments <- guide_assignments %>% 
  add_count(panel, cell, name = "guides_per_cell") %>% 
  filter(guides_per_cell == 1)

# calculate number of cells per target
cells_per_target <- guide_assignments %>% 
  filter(grepl("^chr.+", grna_target)) %>% 
  count(panel, grna_target, name = "cells")

# calculate mean and media cells per target
avg_cells_per_target <- cells_per_target %>% 
  summarize(mean_cells_per_target = mean(cells_per_target$cells),
            median_cells_per_target = median(cells_per_target$cells))

# write to output file
write("# Guide assignment summary statistics", file = snakemake@output[[1]])
write_tsv(summary_stats_across_panels, file = snakemake@output[[1]], col_names = TRUE,
          append = TRUE)
write("\n# Average cells per target (only cells with one assigned guide)",
      file = snakemake@output[[1]], append = TRUE)
write_tsv(avg_cells_per_target, file = snakemake@output[[1]], col_names = TRUE, append = TRUE)
