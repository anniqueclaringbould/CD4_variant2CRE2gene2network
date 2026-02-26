## Get gene constraints for all genes

library(data.table)
library(dplyr)

cs_file <- "/g/steinmetz/project/otar/interpretation/gene_constraints/data/gnomad.v4.1.constraint_metrics.tsv"
cs <- fread(cs_file)

# extract 'Probability of loss-of-function intolerance' for each gene
cs <- cs %>%
  filter(canonical == TRUE) %>%
  select(gene_id, lof.pLI) %>%
  distinct()

# save to output file
fwrite(cs, file = here::here("resources/gene_constraints.csv.gz"))
