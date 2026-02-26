## Create yaml config file containing all enhancer screen guide assignment and sceptre input files

library(yaml)
library(here)

# get all directories containig sceptre data per panel
base_dir <- "/g/steinmetz/moonen/Rstudio/screen_analyses_DM_with_promoterC/results"
sceptre_dirs <- file.path(base_dir, paste0("Panel", 1:16))

# guide assignment matrices
guide_assignments <- file.path(sceptre_dirs, "all_sceptre_output/grna_assignment_matrix.rds")
names(guide_assignments) <- sub(".+(Panel[[:digit:]]+).+", "\\1", guide_assignments)

# get all sceptre objects
sceptre_objects <-  file.path(sceptre_dirs, "sceptre_obj.rds")
names(sceptre_objects) <- sub(".+(Panel[[:digit:]]+).+", "\\1", sceptre_objects)

# get all sceptre objects with full results
finished_sceptre_objects <- list.files(sceptre_dirs, "sceptre_obj_finished", full.names = TRUE, recursive = TRUE)
names(finished_sceptre_objects) <- sub(".+(Panel[[:digit:]]+).+", "\\1", finished_sceptre_objects)

# convert to lists to write to yaml files
guide_assignments <- list(guide_assignments = as.list(guide_assignments))
sceptre_objects <- list(sceptre_objects = as.list(sceptre_objects))
finished_sceptre_objects <- list(finished_sceptre_objects = as.list(finished_sceptre_objects))

# write to yaml file
write(as.yaml(guide_assignments), file = here("config/enhancer_screen_sceptre_files.yml"))
write(as.yaml(sceptre_objects), file = here("config/enhancer_screen_sceptre_files.yml"),
      append = TRUE)
write(as.yaml(finished_sceptre_objects), file = here("config/enhancer_screen_sceptre_files.yml"),
      append = TRUE)
