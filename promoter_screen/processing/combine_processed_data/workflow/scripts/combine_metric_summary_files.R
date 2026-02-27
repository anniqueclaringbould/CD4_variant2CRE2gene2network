## Process and combine QC metrics for each sample from metrics summary files

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

## Define functions --------------------------------------------------------------------------------

# load tables in one metric summary file
process_metric_file <- function(x, tables) {
  
  message("Processing file: ", x)
  
  # load all lines in metric summary file
  lines <- readLines(x)
  
  # load all tables in file
  metrics <- lapply(structure(tables, names = tables), FUN = extract_table_file, lines = lines)
  
  return(metrics)
  
}

# helper function to extract one table from the lines of a metric summary file
extract_table_file <- function(lines, header, show_col_types = FALSE) {
  
  # get lines where table starts based on where the header is found
  table_start <- which(str_detect(lines, paste0("^#", header, "#"))) + 1
  lines <- lines[table_start:length(lines)]
  
  # get end of table based on the first empty line or last line
  lines <- c(lines, "")
  table_end <- which(lines == "")[[1]] - 1
  lines <- lines[1:table_end]
  
  # convert lines to table
  table <- read_csv(paste(lines, collapse = "\n"), show_col_types = show_col_types)
  
  return(table)
  
}

# function to extract and combine specified table from list of tables per sample
combine_table_list <- function(metrics_list, table) {
  
  # extract specified table for each sample and combine into one big table
  combined_table <- metrics_list %>% 
    lapply(FUN = function(x) x[[table]] ) %>% 
    bind_rows(.id = "sample")
  
  return(combined_table)

}

## Load and merge all metric summary files ---------------------------------------------------------

# get all input metric summary files
input_files <- unlist(snakemake@input)
names(input_files) <- basename(dirname(input_files))

# all metric tables to extract and combine
metric_table_names <- structure(snakemake@params$metrics, names = snakemake@params$metrics)

# load and process all specified metric tables in these files
metrics <- lapply(input_files, FUN = process_metric_file, tables = metric_table_names)

# extract and combine all specified metric tables from all samples
metrics <- lapply(metric_table_names, FUN = combine_table_list, metrics_list = metrics)

# save metric tables to output file
saveRDS(metrics, file = snakemake@output[[1]])
