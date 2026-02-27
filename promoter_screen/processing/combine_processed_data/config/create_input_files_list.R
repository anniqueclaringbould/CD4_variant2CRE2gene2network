## Create list of samples for input

library(yaml)

# get all seurat and metric files in input directory
indir <- "/g/steinmetz/project/otar_2063/promoter_screen/processing/screen_data/results"
seurat_files <- list.files(indir, pattern = "Seurat.rds", recursive = TRUE, full.names = TRUE)
metric_files <- list.files(indir, pattern = "Metrics_Summary", recursive = TRUE, full.names = TRUE)

# make simple list of metric files per sample
names(metric_files) <- basename(dirname(metric_files))

# make list of seurat files per donor and chip
donors <- sub("(.)_.+", "\\1", basename(dirname(seurat_files)))
seurat_files <- split(seurat_files, f = donors)

# combine into one list and write to yaml file that can be used like a config file
files <- list(metric_files = metric_files, seurat_files = seurat_files)
write_yaml(files, file = "config/input_files.yml")
