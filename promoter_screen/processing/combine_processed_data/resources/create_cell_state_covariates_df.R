
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(here)
})

# load new cell states covariate table
covars_file <- "/g/steinmetz/project/otar_2063/promoter_screen/interpretation/cell_scoring/output/Representative_AUC_scores_per_cell.tsv"
covars <- read_tsv(covars_file, col_types = cols(cell = col_character(), .default = col_double()))

# reformat for sceptre pipeline
colnames(covars) <- gsub("\\.", "_", colnames(covars))
covars <- covars %>% 
  rename(cell_barcode = cell) %>% 
  as.data.frame()

# write to output file
saveRDS(covars, file = here("resources/cell_state_covariates.rds"))
