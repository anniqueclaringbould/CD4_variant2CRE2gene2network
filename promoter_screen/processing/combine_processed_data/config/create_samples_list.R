## Create list of samples for input

library(yaml)

# get all samples from input directory
indir <- "/g/steinmetz/project/otar_2063/promoter_screen/processing/screen_data/results"
samples <- list.files(indir)

# write to yaml file that can be used like a config file
write_yaml(list("samples" = samples), file = "config/samples.yml")
