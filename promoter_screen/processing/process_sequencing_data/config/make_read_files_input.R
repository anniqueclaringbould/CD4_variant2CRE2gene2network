#####################################################################
## Script to make file with files to process in snakemake pipeline ##
#####################################################################

#### SETUP ####
#libraries
library(data.table)
library(tidyverse)

# Find files in sequencing directories -------------------------------------------------------------

#Find all directories with reads
dirs_embl <- list.dirs(path = "/g/steinmetz/STOCKS/Data/Assay/sequencing/2025", recursive = FALSE)
dirs_embl <- dirs_embl[str_detect(dirs_embl, "Moonen")]  # keep only Dewi's files
dirs_embl <- dirs_embl[!str_detect(dirs_embl, "2025-03")]  # remove the folders with Dewi's files that were not part of this sequencing run
dirs_downtown <- list.dirs(path = "/g/steinmetz/STOCKS/User_Dropboxes/moonen/demultiplexing_done", recursive = FALSE)

#Add the next part of the address
dirs_embl_full <- paste0(dirs_embl, "/Samples/Moonen")
dirs_downtown_full <- paste0(dirs_downtown, "/Demultiplexed/Samples/Moonen")

dirs <- c(dirs_embl_full, dirs_downtown_full)

# Directories with different path patterns that needs to be changed
different_path_dirs <- c(
  "/g/steinmetz/STOCKS/Data/Assay/sequencing/2025/2025-04-12_AV233002_Moonen_2430435701/Samples/Moonen",
  "/g/steinmetz/STOCKS/Data/Assay/sequencing/2025/2025-04-13_AV233002_Moonen_2421466698/Samples/Moonen",
  "/g/steinmetz/STOCKS/Data/Assay/sequencing/2025/2025-04-13_AV233002_Moonen_2421493191/Samples/Moonen"
)

# Manually create new paths for these directories
manual_subdirs <- c("20250412_AV233002_Moonen-25s000392-9", "20250413_AV233002_Moonen-25s000392-12",
                    "20250413_AV233002_Moonen-25s000392-11")
manual_paths <- paste(dirname(dirname(different_path_dirs)), manual_subdirs,
                      "Demultiplexed/Samples/Moonen", sep = "/")

# Replace paths in 'dirs'
dirs <- c(setdiff(dirs, different_path_dirs), manual_paths)

#Get all files
file <- unlist(lapply(dirs, function(d) {
  list.files(path = d, pattern = "\\.fastq\\.gz$", recursive = TRUE, full.names = TRUE)
}))

#Extract sample ID
read_files <- data.frame(file) %>%
  mutate(sample = sub(".*/Moonen/([^/]{4}).*", "\\1", file),
         sample = sub("-", "_", sample)) %>%
  select(sample, file) %>%
  arrange(sample)

# Add files for sample B2_1 from pilot sequencing run ----------------------------------------------

# Directory containing trimmed pilot fastq files
pilot_dir <- "/g/steinmetz/project/otar_2063/promoter_screen/processing/trimmed_fastqs_pilot"

# Get all trimmed fastq files without dump files from trimming
pilot_files <- list.files(path = pilot_dir, pattern = "\\.fastq\\.gz$", recursive = TRUE, full.names = TRUE)
pilot_files <- grep(pilot_files, pattern = "unpaired", invert = TRUE, value = TRUE)

# add to reads table
read_files <- bind_rows(read_files, tibble(sample = "B2_1", file = pilot_files))

# Count the number of files per sample to check if everything is ok
files_per_sample <- count(read_files, sample, name = "files")

# Write output to files
fwrite(read_files, "/g/steinmetz/project/otar_2063/promoter_screen/processing/screen_data/config/read_files_generated.tsv", col.names = TRUE, row.names = FALSE, sep = "\t")
fwrite(files_per_sample, "/g/steinmetz/project/otar_2063/promoter_screen/processing/screen_data/config/read_files_per_sample.csv")
