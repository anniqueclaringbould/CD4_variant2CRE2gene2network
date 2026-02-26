## Combine enhancer screen guide assignments into one file

# save.image("RDA/combine_guide_assignments.rda")
# stop()

suppressPackageStartupMessages(
  library(data.table)
)

# all guide assignment files
guide_assignment_files <- unlist(snakemake@input$guide_assignments)
names(guide_assignment_files) <- sub(".*(Panel[[:digit:]]+).*", "\\1", guide_assignment_files)

# load all guide assignment files and combine into one table
guide_assignments <- lapply(guide_assignment_files, FUN = fread)
guide_assignments <- rbindlist(guide_assignments, idcol = "panel")

# all summary statistics files
summary_stats_files <- unlist(snakemake@input$summary_stats)
names(summary_stats_files) <- sub(".*(Panel[[:digit:]]+).*", "\\1", summary_stats_files)

# load all summary statistics files and combine into one table
summary_stats <- lapply(summary_stats_files, FUN = fread)
summary_stats <- rbindlist(summary_stats, idcol = "panel")

# all grna targets table files
grna_targets_files <- unlist(snakemake@input$grna_targets)
names(grna_targets_files) <- sub(".*(Panel[[:digit:]]+).*", "\\1", grna_targets_files)

# load all guide assignment files and combine into one table
grna_targets <- lapply(grna_targets_files, FUN = fread)
grna_targets <- rbindlist(grna_targets, idcol = "panel")

# save to output file
fwrite(guide_assignments, file = snakemake@output$guide_assignments, quote = FALSE, sep = "\t", na = "NA")
fwrite(summary_stats, file = snakemake@output$summary_stats, quote = FALSE, sep = "\t", na = "NA")
fwrite(grna_targets, file = snakemake@output$grna_targets, quote = FALSE, sep = "\t", na = "NA")
