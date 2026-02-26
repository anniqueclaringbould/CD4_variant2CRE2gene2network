## Create yaml file containing all inputs for the BD Rhapsody pipeline for one sample

# required packages
suppressPackageStartupMessages({
  library(yaml)
  library(R.utils)
})

# Define functions ---------------------------------------------------------------------------------

# get default snakemake workflow directory
get_workflow_dir <- function() {
  dirname(dirname(snakemake@scriptdir))
}

# helper function to convert files to absolute paths if they are relative paths
absolute_path <- function(file, workflow_dir = get_workflow_dir()) {
  as.character(getAbsolutePath(file, workDirectory = workflow_dir))
}

# function to create entry for one reads file
make_reads_entry <- function(reads_file) {
  list(
    class = "File",
    location = reads_file
  ) 
}

# Create input yaml file ---------------------------------------------------------------------------

# create list for all input reads files for the given sample
read_files <- list(
  Reads = lapply(absolute_path(snakemake@input$reads), FUN = make_reads_entry)
)

# create list for reference
reference <- list(
  Reference_Archive = list(
    class = "File",
    location = absolute_path(snakemake@input$reference)
  )
)

# create entry for optional cell numbers if provided
if (!is.na(snakemake@params$cell_numbers)) {
  cells <- list(
    Expected_Cell_Count = as.integer(snakemake@params$cell_numbers)
  )
} else {
  cells <- NULL
}

# create list for other pipeline arguments
arguments <- list(
  Run_Name = snakemake@wildcards$sample,
  Generate_Bam = verbatim_logical(snakemake@params$generate_bam),
  Maximum_Threads = as.integer(snakemake@params$threads)
)

# combine into list specifying all inputs and arguments
all_input <- c(
  "cwl:tool" = "rhapsody",
  read_files,
  reference,
  cells,
  arguments
)

# conver to yamls string without trailing newline character
input_yaml <- sub("\n$", "", as.yaml(all_input))

# write to output yaml file
write("#!/usr/bin/env cwl-runner", file = snakemake@output[[1]])
write(input_yaml, file = snakemake@output[[1]], append = TRUE)
