## Extract significant promoter screen hits from DESeq2 pseudobulked differential expression results

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(here)
})

# list of promoter screen output files
prom_results_dir <- "path/to/deseq2/output"
prom_results_files <- list.files(prom_results_dir, full.names = TRUE)

# function to load one results file and extract significant hits
load_file <- function(file, fdr_threshold = 0.1) {
  x <- file %>% 
    fread() %>% 
    filter(adj_p_value < fdr_threshold)
  return(x)
}

# load all files
prom_results <- lapply(prom_results_files, FUN = load_file, fdr_threshold = 0.1)
prom_results <- rbindlist(prom_results)

# write to output file
fwrite(prom_results, file = "promoter_screen_hits.tsv.gz", na = "NA")
