## Create on-disc sceptre object from counts data in matrix (.mtx) files, which emulate data
## processed with cell ranger

# save.image("create_sceptre_object.rda")
# stop()

# opening log file to collect all messages, warnings and errors
log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type = "message")

# required packages
suppressPackageStartupMessages({
  library(sceptre)
})

# Prepare input ------------------------------------------------------------------------------------

message("Loading input files...")

# load gRNA targets table
grna_targets <- read.table(snakemake@input$grna_targets, header = TRUE, sep = "\t")

# load file with additional, pre-computed cell covariates
cell_covariates <- readRDS(snakemake@input$covariates)

# load positive control pairs if provided
if (!is.null(snakemake@input$pos_ctrl_pairs)) {
  pos_ctrl_pairs <- read.table(snakemake@input$pos_ctrl_pairs, header = TRUE, sep = "\t")
} else {
  pos_ctrl_pairs <- data.frame(grna_target = character(0), response_id = character(0))
}

# get all directories containing count data
count_dirs <- dirname(snakemake@input$mtx)

# load all cell barcodes in specified input files
cell_barcodes <- unlist(lapply(file.path(count_dirs, "barcodes.tsv.gz"), FUN = readLines))

# filter cell covariates table for used cell barcodes only
cell_covariates <- cell_covariates[rownames(cell_covariates) %in% cell_barcodes, ]

# subset cell covariates table to covariates to include if provided in parameters
if (!is.null(snakemake@params$covars)) {
  cell_covariates <- cell_covariates[, snakemake@params$covars] 
}

# Create sceptre object ----------------------------------------------------------------------------

message("Creating sceptre object:")

# create sceptre object using on-disc data structure
sceptre_object <- import_data_from_cellranger(
  directories = count_dirs,
  moi = snakemake@params$moi,
  grna_target_data_frame = grna_targets,
  extra_covariates = cell_covariates,  
  use_ondisc = snakemake@params$on_disc,
  directory_to_write = dirname(snakemake@output[[1]])
)

# TODO: FIX THIS
# remove batch covariate if not in formula, since else sceptre doesn't work properly
if (!grepl("batch", snakemake@params$formula)) {
  sceptre_object@covariate_data_frame <- sceptre_object@covariate_data_frame[, -ncol(sceptre_object@covariate_data_frame)]
  sceptre_object@covariate_names <- setdiff(sceptre_object@covariate_names, "batch")
}

message("Setting analysis parameters...")

# set analysis parameters
sceptre_object <- set_analysis_parameters(
  sceptre_object = sceptre_object,
  positive_control_pairs = pos_ctrl_pairs,
  side = snakemake@params$side,
  formula_object = as.formula(snakemake@params$formula),
  resampling_mechanism = snakemake@params$resampling_mechanism
)

# Save sceptre object to output file ---------------------------------------------------------------

message("Writing sceptre object to file...")

# save sceptre object to file
if (snakemake@params$on_disc == TRUE) {
  write_ondisc_backed_sceptre_object(
    sceptre_object = sceptre_object,
    directory_to_write = dirname(snakemake@output[[1]])
  )
} else {
  saveRDS(sceptre_object, file = snakemake@output[[1]])
}

message("All done!")

# close log file connection
sink()
sink(type = "message")
close(log)
