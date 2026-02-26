## Extract trans-acting effects from promoter screen results

# save.image("RDA/extractTransProm.rda")
# stop()

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# load processed sceptre output
sig_trans_pairs <- read_tsv(snakemake@input$results, show_col_types = FALSE)

# reformat for output
sig_trans_pairs <- sig_trans_pairs %>% 
  select(promoter = grna_target, trans_gene = gene_symbol) %>% 
  arrange(promoter, trans_gene)

# save to output file
write_csv(sig_trans_pairs, file = snakemake@output[[1]])
