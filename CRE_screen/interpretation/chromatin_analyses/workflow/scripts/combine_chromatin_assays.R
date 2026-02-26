## Combine all enhancer chromatin assay signals and add to results table

# save.image("RDA/combine_chromatin_assays.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
})

## Define functions --------------------------------------------------------------------------------

# load and combine assay count files into one table
combine_assay_counts <- function(assay_files, assay_names, tag = NULL) {
  
  # load all assay count files and convert to one table
  assay_counts <- assay_files %>% 
    lapply(FUN = read_tsv, show_col_types = FALSE) %>% 
    bind_rows(.id = "file")
  
  # add chromatin assay name to table
  assay_counts <- left_join(assay_counts, assay_names, by = "file")
  
  # normalize signal by element size to calculate average signal across element (for FC measurements)
  assay_counts <- mutate(assay_counts, length = end - start, signal_avg = signal_count / length)
  
  # convert to wide format
  assay_counts <- assay_counts %>% 
    select(-c(signal_total, file, length)) %>% 
    pivot_wider(values_from = c(signal_count, signal_avg), names_from = assay)
  
  # add tag to chromatin signal columns identifying these as signal at enhancers
  assay_counts <- rename_with(assay_counts, ~ paste0(tag, .x), starts_with("signal_"))
  
  return(assay_counts)
  
}

## Combine and add all enhancer assays -------------------------------------------------------------

# load CRISPR results
results <- read_csv(snakemake@input$elements, show_col_types = FALSE)

# get table with chromatin assay names for all assay files from ENCODE 
assay_names <- enframe(unlist(snakemake@config$chromatin_assays), name = "assay", value = "file")

# files containing assay counts
enh_assay_files <- snakemake@input$enh_assays
names(enh_assay_files) <- sub(".+enhancers\\.(.+)\\.tsv", "\\1", enh_assay_files)

# load and combine all enhancer chromatin assay files
enh_assay_counts <- combine_assay_counts(enh_assay_files, assay_names = assay_names, tag = "enh_")

# merge with CRISPR results
results <- enh_assay_counts %>% 
  select(-c(chr, start, end)) %>% 
  left_join(results, ., by = c("grna_target" = "name"))

## Combine and add all promoter assays -------------------------------------------------------------

# files containing assay counts
prom_assay_files <- snakemake@input$prom_assays
names(prom_assay_files) <- sub(".+promoters\\.(.+)\\.tsv", "\\1", prom_assay_files)

# load and combine all promoter chromatin assay files
prom_assay_counts <- combine_assay_counts(prom_assay_files, assay_names = assay_names,
                                          tag = "tss_")

# merge with CRISPR results
results <- prom_assay_counts %>% 
  select(-c(chr, start, end)) %>% 
  left_join(results, ., by = c("transcript_tss" = "name"))

# save to output file
write_csv(results, snakemake@output[[1]])
