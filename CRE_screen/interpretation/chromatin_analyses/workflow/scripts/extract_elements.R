## Extract enhancer or promoter elements from E-G results and save to bed file

# save.image("RDA/extract_elements.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

# load all tested E-G pairs
pairs <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

# extract unique enhancer or promoter elements and convert to table in bed format
if (snakemake@params$elements == "enhancers") {
  
  elements <- pairs %>% 
    filter(!is.na(pert_chr)) %>% 
    select(pert_chr, pert_start, pert_end, grna_target) %>% 
    distinct()
  
} else if(snakemake@params$elements == "promoters") {
  
  # get how much each tss needs to be extended up- and downstream
  extend <- snakemake@params$extend_tss
  if (is.null(extend)) {
    stop("'extend_tss' parameter required for extracting promoters", call. = FALSE)
  }
  
  # extract gene TSSs and extend to create promoter regions
  elements <- pairs %>% 
    filter(!is.na(gene_tss)) %>% 
    mutate(tss_start = gene_tss - extend,
           tss_end = gene_tss + extend) %>% 
    select(gene_chr, tss_start, tss_end, transcript_tss) %>% 
    distinct()
  
} else {
  stop("Invalid 'elements' parameter", call. = FALSE)
}

# save to output file
write_tsv(elements, file = snakemake@output[[1]], col_names = FALSE)
