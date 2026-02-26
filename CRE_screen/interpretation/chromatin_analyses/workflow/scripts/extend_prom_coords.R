## Extend manual coordinates of control promoters by given amount for chromatin quantifications

# save.image("RDA/extend_prom_coords.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# load manually annotated perturbation coordinates for promoter controls
pert_coords <- read_tsv(snakemake@input[[1]], show_col_types = FALSE)

# extend coordinates of promoter targeting controls by specified amount around TSS
extend <- snakemake@params$extend_tss
pert_coords <- pert_coords %>% 
  mutate(
    pert_start = if_else(ctrl_type == "promoter", true = pert_start - extend, false = pert_start),
    pert_end = if_else(ctrl_type == "promoter", true = pert_end + extend, false = pert_end)
    )

# save extended coordinates to output file
write_tsv(pert_coords, file = snakemake@output[[1]])
