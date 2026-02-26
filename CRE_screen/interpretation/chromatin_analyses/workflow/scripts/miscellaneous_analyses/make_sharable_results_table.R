## Make results table to share with OTAR

# save.image("RDA/make_sharable_results_table.rda")
# stop()

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

# load full enhancer screen results table
etp <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

# for control perturbations which were in multiple panels, select the most significant one
etp_ctrl <- etp %>% 
  filter(ctrl_perturbation == TRUE) %>% 
  group_by(response_id, grna_target) %>% 
  slice_min(p_value, n = 1, with_ties = FALSE)

# add back to results table with discovery pairs
etp <- etp %>% 
  filter(ctrl_perturbation == FALSE) %>% 
  bind_rows(etp_ctrl) %>% 
  arrange(desc(p_value))

# get all feature columns to remove
chrom_assay_cols <- grep(colnames(etp), pattern = "enh_signal|tss_signal", value = TRUE)
hic_cols <- "hic_int_freq"
gene_features <- c("gene_constraints", "gene_uniq_expr", "gene_subcellularLocations",
                   "gene_drugIds", "gene_maxClinicalPhase", "gene_references", "gene_diseaseIds",
                   "gene_diseaseLabels", "gene_therapeuticAreas", "gene_clinicalAvailable")

# remove all these columns from the table
etp <- select(etp, -all_of(c(chrom_assay_cols, hic_cols, gene_features)))

# reorder columns and ensure correct sorting
etp <- etp %>% 
  relocate(grna_target, .before = 1) %>% 
  relocate(gene_tpm, .after = "gene_tss") %>% 
  arrange(p_value, grna_target, response_id)

# save to output file
write_csv(etp, snakemake@output[[1]])
