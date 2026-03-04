
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install necessary Bioconductor packages
BiocManager::install(c("TFBSTools", "JASPAR2022", "Biostrings",
                       "GenomicRanges", "BSgenome", "rtracklayer", "org.Hs.eg.db"))

if (!requireNamespace("motifmatchr", quietly = TRUE))
    BiocManager::install("motifmatchr")


library(circlize)
library(tidyverse)
library(data.table)
library(Seurat)
library(optparse)
library(ggplot2)
library(ggrepel)
library(viridis)
library(gridExtra)
library(Rcpp)
library(scDblFinder)
library(ggrastr)
library(Matrix)
library(tibble)
library(dplyr)
library(RColorBrewer)
library(ggpubr)
library(clusterProfiler)
library(rtracklayer)
library(motifmatchr)
library(TFBSTools)
library(JASPAR2022)
library(Biostrings)
library(GenomicRanges)
library(BSgenome)
library(rtracklayer)
library(org.Hs.eg.db)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(AnnotationDbi)
library(AnnotationHub)
library(readxl)

select <- dplyr::select
rename <- dplyr::rename
slice <- dplyr::slice
map <- purrr::map
mutate <- dplyr::mutate
filter<-dplyr::filter

outdir <- paste0("/g/steinmetz/moonen/Rstudio/screen_analyses_DM/results")

results_df <- read_csv('/g/steinmetz/moonen/Rstudio/screen_analyses_DM/results/results_df_annotated.csv')

tpm_data <- read.csv("/g/steinmetz/project/otar/promoter_screen/output/gene_list_unfiltered.txt", sep = "\t")

# TRRUST <- read_tsv('/g/steinmetz/moonen/Rstudio/screen_analyses_DM/results/TFBS/trrust_rawdata.human.tsv')
# colnames(TRRUST) <- c("TF", "Target", "Regulation", "PMID")
# repressive_tf_symbols <- TRRUST %>% 
#   filter(Regulation == "Repression") %>% 
#   select(TF)%>%
#   pull(TF) %>%
#   unique()
# activator_tf_symbols <- TRRUST %>% 
#   filter(Regulation == "Activation") %>% 
#   select(TF) %>%
#   pull(TF) %>%
#   unique()

file_path <- "/g/steinmetz/moonen/Rstudio/screen_analyses_DM/resources/1-s2.0-S2211124719314391-mmc2.xls"
sheet_names <- excel_sheets(file_path)
print(sheet_names)
# Create empty lists to store the results
activators_list <- list()
repressors_list <- list()

for (i in 2:length(sheet_names)) {
  sheet <- sheet_names[i]
  
  # Read the current sheet
  df <- read_excel(file_path, sheet = sheet, col_names = TRUE)
  print(colnames(df))
  print(head(df))
  print(head(df$classification))
  # Check if the classification column exists in the dataframe
  if ("classification" %in% colnames(df)) {
    # Filter for activators and take only the TF column
    activators <- df[df$classification == "activator", "TF"]
    
    # Filter for repressors and take only the TF column
    repressors <- df[df$classification == "repressor", "TF"]
    
    # Append the TF values to the lists
    activators_list[[sheet]] <- activators
    repressors_list[[sheet]] <- repressors
  } else {
    # Print a message for sheets without the classification column
    print(paste("Skipping sheet:", sheet, "- no classification column"))
  }
}

# Combine the lists into single dataframes
activator_tf_symbols <- do.call(rbind, activators_list)
activator_tf_symbols <- activator_tf_symbols %>% distinct(TF)
repressive_tf_symbols <- do.call(rbind, repressors_list)
repressive_tf_symbols <- repressive_tf_symbols %>% distinct(TF)


overlap <- intersect(activator_tf_symbols, repressive_tf_symbols)
# Remove the overlapping elements from the original dataframes
activator_tf_symbols <- setdiff(activator_tf_symbols, overlap)
repressive_tf_symbols <- setdiff(repressive_tf_symbols, overlap)

print(head(repressive_tf_symbols))
print(head(activator_tf_symbols))
print(length(repressive_tf_symbols))
print(length(activator_tf_symbols))

# Identify rows with genomic coordinates in 'grna_target'
has_coords <- grepl("^chr[0-9XYM]+:[0-9]+-[0-9]+$", results_df$grna_target)

# Create a new data frame 'results_df_coords' containing only entries with genomic coordinates
results_df_coords <- results_df[has_coords, ]

# Reset row names (optional)
rownames(results_df_coords) <- NULL


# Use 'str_match' to extract chromosome, start, and end positions
parsed_coords <- str_match(results_df_coords$grna_target, "^([^:]+):([0-9]+)-([0-9]+)$")

# The 'parsed_coords' matrix has columns:
# 1. Full match
# 2. Chromosome
# 3. Start
# 4. End

# Add new columns to 'results_df_coords' without modifying 'results_df'
results_df_coords$chr <- parsed_coords[,2]
results_df_coords$start <- as.numeric(parsed_coords[,3])
results_df_coords$end <- as.numeric(parsed_coords[,4])

results_df_coords <- results_df_coords %>% distinct(grna_target, .keep_all = TRUE)
# Create a GRanges object 'regions' from 'results_df_coords'
regions <- GRanges(seqnames = results_df_coords$chr,
                   ranges = IRanges(start = results_df_coords$start, end = results_df_coords$end),
                   strand = "*",
                   mcols = results_df_coords)


genome <- BSgenome.Hsapiens.UCSC.hg38
sequences <- getSeq(genome, regions)


# Set options to retrieve vertebrate TFs (species: Homo sapiens, taxon ID 9606)
#pfm_directory <- "/g/steinmetz/moonen/Rstudio/screen_analyses_DM/resources/JASPAR2024_CORE_vertebrates_non-redundant_pfms_jaspar/"
pfm_directory <- "/g/steinmetz/moonen/Rstudio/screen_analyses_DM/resources/pfm/"
pfm_files <- list.files(pfm_directory, pattern = "\\.jaspar$", full.names = TRUE)


# Read each PFM file into a list
pfm_list <- lapply(pfm_files, function(file) {
    pfms <- readJASPARMatrix(file, matrixClass = "PFM")
    # Check if the result is a PFMatrixList
    if (is(pfms, "PFMatrixList")) {
        # Convert PFMatrixList to a list of PFMatrix objects
        return(as.list(pfms))
    } else {
        # Wrap PFMatrix in a list
        return(list(pfms))
    }
})

pfm_list <- unlist(pfm_list, recursive = FALSE)

# Assign names to the list elements using the 'name' slot of PFMatrix
names(pfm_list) <- sapply(pfm_list, function(pfm) pfm@name)


# go_ids_repressor <- c(
#     "GO:0000119",  # Transcription Repressor Activity
#     "GO:0000124",  # Negative Regulation of Transcription from RNA Polymerase II Promoter
#     "GO:0000125",  # Negative Regulation of Transcription, DNA-Templated
#     "GO:0006357",  # Negative Regulation of Transcription by RNA Polymerase II 
#     "GO:0030529",  # Repressor Activity, RNA Polymerase II Promoter
#     "GO:0106250",  # Negative Regulation of Transcription from RNA Polymerase II Promoter, Initiation
#     "GO:0001217",  # Negative Regulation of Transcription by RNA Polymerase I
#     "GO:0001227",   # Negative Regulation of Transcription by RNA Polymerase III
#     "GO:0016564", #Transcription repressor complex
#     "GO:0000122", #Negative regulation of transcription by RNA polymerase II
#     "GO:0016573" #Histone deacetylation
# )

# go_ids_repressor <- c(
#     "GO:0003714",  # Transcription Corepressor Activity
#     "GO:0016564",  # Transcription Repressor Complex
#     "GO:0016573",  # Histone Deacetylation
#     "GO:0071549",  # Chromatin Silencing
#     "GO:0000122",  # Negative Regulation of Transcription by RNA Polymerase II
#     "GO:0045892",  # Negative Regulation of Transcription, DNA-templated
#     "GO:0032968",  # Positive Regulation of Histone Deacetylation
#     "GO:0000121",  # Negative Regulation of Transcription from RNA Polymerase II Promoter, Initiation
#     "GO:0000125",  # Negative Regulation of Transcription, DNA-templated
#     "GO:0003713",  # Transcription Repressor Activity, RNA Polymerase II Core Promoter Proximal Region
#     "GO:0090071"   # Negative Regulation of Histone H3-K4 Methylation
# )

# go_ids_repressor<- c(
#     "GO:0001227",  # DNA-binding transcription repressor activity, RNA polymerase II-specific
#     "GO:0003713",  # Transcription repressor activity, RNA polymerase II core promoter proximal region
#     "GO:0017053",  # Transcription repressor binding
#     "GO:0000977",  # RNA polymerase II core promoter sequence-specific DNA binding, negative regulation
#     "GO:0000122",  # Negative regulation of transcription by RNA polymerase II
#     "GO:0003714",  # Transcription corepressor activity
#     "GO:0016564",  # Transcription repressor complex
#     "GO:0030529",  # Repressor activity, RNA polymerase II promoter
#     "GO:0003712",   # Negative regulation of transcription core promoter proximal region sequence-specific DNA binding
#     "GO:0071549",  # Chromatin silencing
#     "GO:0016573",  # Histone deacetylation
#     "GO:0090071",  # Negative regulation of histone H3-K4 methylation
#     "GO:0000121"  # Negative regulation of transcription from RNA polymerase II promoter, initiation
# )

# go_ids_repressor<- TRRUST %>% filter(regulation=="Repression") %>% select(Gene)

# # Retrieve all genes associated with the GO terms
# repressive_genes <- AnnotationDbi::select(org.Hs.eg.db,
#                                           keys = go_ids_repressor,
#                                           columns = c("SYMBOL"),
#                                           keytype = "GOALL")

# # Get unique gene symbols
# repressive_tf_symbols <- unique(repressive_genes$SYMBOL)


# go_ids_activator <- c(
#     "GO:0001223",  # Transcription Coactivator Activity
#     "GO:0010468",  # RNA Polymerase II Transcription Coactivator Activity
#     "GO:0001227",  # DNA-binding Transcription Activator Activity
#     "GO:0001228",  # Positive Regulation of Transcription from RNA Polymerase II Promoter
#     "GO:0045944",  # Positive Regulation of Transcription by RNA Polymerase II
#     "GO:0045893",  # Positive Regulation of Transcription, DNA-Templated
#     "GO:0045945"   # Positive Regulation of Transcription from RNA Polymerase II Promoter
# )

# go_ids_activator<- c(
#     "GO:0001228",  # DNA-binding transcription activator activity, RNA polymerase II-specific
#     "GO:0000978",  # RNA polymerase II cis-regulatory region sequence-specific DNA binding
#     "GO:0003713",  # Transcription activator activity, RNA polymerase II core promoter proximal region
#     "GO:0001077",  # RNA polymerase II core promoter proximal region sequence-specific DNA binding, positive regulation
#     "GO:0006355",  # Regulation of transcription, DNA-templated
#     "GO:0001227",  # Transcription coactivator activity
#     "GO:0035328",  # Enhancer sequence-specific DNA binding
#     "GO:0001225"   # DNA-binding transcription activator activity, enhancer binding
# )


# activator_genes <- AnnotationDbi::select(org.Hs.eg.db,
#                                           keys = go_ids_activator,
#                                           columns = c("SYMBOL"),
#                                           keytype = "GOALL")

# # Get unique gene symbols
# activator_tf_symbols <- unique(activator_genes$SYMBOL)




repressive_tf_symbols <- repressive_tf_symbols %>%
  # Identify genes in tpm_data with TPM >= 3
  intersect(tpm_data %>% filter(tpm >= 1) %>% pull(Gene)) %>%
  # Ensure uniqueness (optional if activator_tf_symbols is already unique)
  unique()

print(head(repressive_tf_symbols))


# Inspect the list
length(repressive_tf_symbols)

# Filter activator_tf_symbols to keep only those with TPM >= 3
activator_tf_symbols <- activator_tf_symbols %>%
  # Identify genes in tpm_data with TPM >= 3
  intersect(tpm_data %>% filter(tpm >= 1) %>% pull(Gene)) %>%
  # Ensure uniqueness (optional if activator_tf_symbols is already unique)
  unique()

length(activator_tf_symbols)



pfm_tfs <- names(pfm_list)

# Match PFMs with repressive TF symbols
repressive_pfm_list <- pfm_list[pfm_tfs %in% repressive_tf_symbols]

# Check how many PFMs you have after filtering
length(repressive_pfm_list)

# List the TFs included after filtering
names(repressive_pfm_list)



# Convert both to upper case
pfm_tfs_upper <- toupper(pfm_tfs)
repressive_tf_symbols_upper <- toupper(repressive_tf_symbols)

# Filter PFMs
repressive_pfm_list <- pfm_list[pfm_tfs_upper %in% repressive_tf_symbols_upper]

# Convert PFMs to PWMs
pwm_list <- lapply(repressive_pfm_list, toPWM)

# Step 1: Extract the first PWMatrix object
pwm_first <- pwm_list[[1]]


# Step 3: Create a PWMatrixList from the entire pwm_list
# Method 1: Using do.call
pwm_set <- do.call(PWMatrixList, pwm_list)


# Verify the PWMatrixList
print(class(pwm_set))
# [1] "PWMatrixList"
print(length(pwm_set))


# Set names
#names(pwm_set) <- names(pwm_list)

# Initialize an empty list to store results
tfbs_results <- list()
# Extract 'grna_target' values
sequence_ids <- results_df_coords$grna_target 


# Ensure that the length of 'sequence_ids' matches the length of 'sequences'
length(sequence_ids)  # Should be equal to length of sequences
length(sequences)

# Assign the 'grna_target' values as names to your sequences
names(sequences) <- sequence_ids

tfbs_results <- list()
# Loop over each sequence
for (i in seq_along(sequences)) {
    seq <- sequences[[i]]             # Get the sequence
    seq_name <- names(sequences)[i]   # Get the sequence name (grna_target)

    # Search for TFBS in the sequence
    hits <- searchSeq(
        pwm_set,                      # Your PWMatrixList object
        seq,                          # The current sequence
        seqname = seq_name,           # The name of the sequence
        min.score = "95%",            # Threshold for motif matching
        strand = "*"                  # Search both strands
    )

    # Store the results in the list using seq_name as the key
    tfbs_results[[seq_name]] <- hits
}

saveRDS(tfbs_results, paste0(outdir, "/TFBS/tfbs_results_repressive_tpm5.RDS"))

# Match PFMs with activator TF symbols
activator_pfm_list <- pfm_list[pfm_tfs %in% activator_tf_symbols]

# Check how many PFMs you have after filtering
length(activator_pfm_list)

# List the TFs included after filtering
names(activator_pfm_list)

activator_tf_symbols_upper <- toupper(activator_tf_symbols)

# Filter PFMs
activator_pfm_list <- pfm_list[pfm_tfs_upper %in% activator_tf_symbols_upper]

# Convert PFMs to PWMs
pwm_list_act <- lapply(activator_pfm_list, toPWM)

# Step 1: Extract the first PWMatrix object
pwm_first_act <- pwm_list_act[[1]]

# Step 3: Create a PWMatrixList from the entire pwm_list
# Method 1: Using do.call
pwm_set_act <- do.call(PWMatrixList, pwm_list_act)


# Verify the PWMatrixList
print(class(pwm_set_act))
# [1] "PWMatrixList"
print(length(pwm_set_act))

# Set names
#names(pwm_set) <- names(pwm_list)

# Initialize an empty list to store results
tfbs_results_act <- list()


# Loop over each sequence
for (i in seq_along(sequences)) {
    seq <- sequences[[i]]             # Get the sequence
    seq_name <- names(sequences)[i]   # Get the sequence name (grna_target)

    # Search for TFBS in the sequence
    hits_act <- searchSeq(
        pwm_set_act,                      # Your PWMatrixList object
        seq,                          # The current sequence
        seqname = seq_name,           # The name of the sequence
        min.score = "95%",            # Threshold for motif matching
        strand = "*"                  # Search both strands
    )

    # Store the results in the list using seq_name as the key
    tfbs_results_act[[seq_name]] <- hits_act
}
saveRDS(tfbs_results_act, paste0(outdir, "/TFBS/tfbs_results_activator_tpm1.RDS"))

# Initialize a vector to store counts
repressive_tfbs_counts <- sapply(tfbs_results, function(hits_list) {
    # Sum the number of hits for all TFs in the sequence
    sum(sapply(hits_list, length))
})

# Assign names to the counts vector
names(repressive_tfbs_counts) <- names(tfbs_results)



# Create a data frame from the counts
tfbs_counts_df_rep <- data.frame(
    grna_target = names(repressive_tfbs_counts),
    repressive_tfbs_count = repressive_tfbs_counts,
    stringsAsFactors = FALSE
)

# Ensure 'grna_target' is of type character in both data frames
results_df$grna_target <- as.character(results_df$grna_target)

# Merge the counts into results_df
results_df <- merge(
    results_df,
    tfbs_counts_df_rep,
    by = "grna_target",
    all.x = TRUE
)


# Initialize a vector to store counts
activator_tfbs_counts <- sapply(tfbs_results_act, function(hits_list) {
    # Sum the number of hits for all TFs in the sequence
    sum(sapply(hits_list, length))
})

# Assign names to the counts vector
names(activator_tfbs_counts) <- names(tfbs_results_act)



# Create a data frame from the counts
tfbs_counts_df_act <- data.frame(
    grna_target = names(activator_tfbs_counts),
    activator_tfbs_count = activator_tfbs_counts,
    stringsAsFactors = FALSE
)

# Merge the counts into results_df
results_df <- merge(
    results_df,
    tfbs_counts_df_act,
    by = "grna_target",
    all.x = TRUE
)

write.csv(results_df, paste0(outdir, "/TFBS/results_df_TFBS_tpm1.csv"), row.names = FALSE)

