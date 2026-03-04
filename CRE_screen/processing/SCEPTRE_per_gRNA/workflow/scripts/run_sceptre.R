## Analyze one OTAR 10x lane using Sceptre

# save.image("run_sceptre.rda")
# stop()

.libPaths(c("/home/moonen/R/x86_64-pc-linux-gnu-library/4.2", .libPaths()))
cat("Using libs:\n"); print(.libPaths())

# required packages
library(ggplot2)
library(readr)
library(dplyr)
library(sceptre, lib.loc="/home/moonen/R/4.2.2-foss-2022b")

# set seed for any random processes
set.seed(20240708)

# load sceptre object
sceptre_obj <- readRDS(snakemake@input$sceptre)

# load gene annotations table
genes <- read.table(snakemake@input$genes, header = TRUE, sep = "\t")

## Define analysis parameters ----------------------------------------------------------------------

# define positive control E-G pairs (promoter controls)
prom_ctrl_pairs <- construct_positive_control_pairs(sceptre_obj)

# define discovery E-G pairs based on provided maximum distance parameter
if (!is.null(snakemake@params$max_dist)) {
  
  # construct cis-pairs within specified distance
  discovery_pairs <- construct_cis_pairs(sceptre_obj,
                                         positive_control_pairs = prom_ctrl_pairs, 
                                         distance_threshold = 1e6,
                                         response_position_data_frame = genes)

} else {
  
  # construct pairs for all perturbation - gene combinations
  discovery_pairs <- construct_trans_pairs(sceptre_obj,
                                           positive_control_pairs = prom_ctrl_pairs,
                                           pairs_to_exclude = "none")
  
}

# remove any duplicated rows due to genes having multiple annotated TSSs
discovery_pairs <- unique(discovery_pairs)

# set analysis parameters
sceptre_obj <- set_analysis_parameters(
  sceptre_object = sceptre_obj,
  discovery_pairs = discovery_pairs,
  positive_control_pairs = prom_ctrl_pairs,
  grna_integration_strategy = "singleton"
)

## Assign gRNAs to cells ---------------------------------------------------------------------------

# assign guides to cell using mixture model approach
sceptre_obj <- assign_grnas(sceptre_obj, method = "mixture", parallel = snakemake@params$parallel)

# plot guide assignment summary statistics plot and visualize 4 randomly picked guides
plot(sceptre_obj, n_grnas_to_plot = 4)
ggsave(file = snakemake@output$guides_plot, width = 9, height = 6)

## Run QC analyses ---------------------------------------------------------------------------------

# plot covariate distributions
plot_covariates(sceptre_obj)
ggsave(file = snakemake@output$covars_plot, width = 7, height = 3.5)

# perform and plot QC analyses
sceptre_obj <- run_qc(sceptre_obj)
plot(sceptre_obj)
ggsave(file = snakemake@output$qc_plot, width = 5, height = 6)

## Run model calibration check ---------------------------------------------------------------------

# run calibration check using non-targeting control guides
sceptre_obj <- run_calibration_check(sceptre_obj, parallel = snakemake@params$parallel)

# plot calibration results
plot(sceptre_obj)
ggsave(file = snakemake@output$calibration_plot, width = 6.5, height = 5)

## Run power check and discovery analysis ----------------------------------------------------------

# run power check using positive control perturbations
sceptre_obj <- run_power_check(sceptre_obj, parallel = snakemake@params$parallel)

# plot power check results
plot(sceptre_obj)
ggsave(file = snakemake@output$power_plot, width = 4, height = 3)

# perform discovery analysis of enhancer targeting perturbations
sceptre_obj <- run_discovery_analysis(sceptre_obj, parallel = snakemake@params$parallel)

# plot discovery analysis results
plot(sceptre_obj)
ggsave(file = snakemake@output$discovery_plot, width = 6, height = 5)

# extract discovery results in table format and save to output file
results <- get_result(sceptre_object = sceptre_obj, analysis = "run_discovery_analysis")
write_tsv(results, file = snakemake@output$results)

# write all output to directory
outdir <- dirname(snakemake@output$results_rds)
write_outputs_to_directory(sceptre_obj, outdir)
