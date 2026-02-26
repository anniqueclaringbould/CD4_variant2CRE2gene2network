## Make input yaml file to generate alignment reference

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

# Create input yaml file ---------------------------------------------------------------------------

# make list of main inputs
input <- list(
  Genome_fasta = list(list(class = "File", location = absolute_path(snakemake@input$genome_fasta))),
  Gtf = list(list(class = "File", location = absolute_path(snakemake@input$genome_gtf)))
)

# add optional extra sequences
if (!is.null(snakemake@input$extra_sequences)) {
  extra_sequences <- list(
    Extra_sequences = list(list(class = "File",
                                location = absolute_path(snakemake@input$extra_sequences)))
  )
  input <- c(input, extra_sequences)
}

# save yaml to file
write_yaml(input, file = snakemake@output[[1]])
