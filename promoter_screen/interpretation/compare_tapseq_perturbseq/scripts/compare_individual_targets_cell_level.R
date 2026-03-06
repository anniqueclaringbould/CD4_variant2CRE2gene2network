####################################################################
## Script to compare TAP-seq and Perturb-seq on single cell level ##
####################################################################

#### SETUP ####
#libraries
library(data.table)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(Seurat)
library(Matrix)

#directories
root <- "/g/steinmetz/project/otar_2063/"
indir.tap <- "/g/steinmetz/moonen/Rstudio/screen_analyses_DM_with_promoterC/results/"
indir.pert <- paste0(root, "promoter_screen/interpretation/seurat_qc/output/combined/")
enrdir <- paste0(root, "promoter_screen/interpretation/enrichment/")
outdir <- paste0(root, "promoter_screen/interpretation/compare_tapseq_perturbseq/output/")


#### READ DATA ####
#TAP-seq
tap_files <- list.files(path = indir.tap, pattern = "merged_seurat_grna_annotated.rds", recursive = T, full.names = T)
tap_list <- lapply(tap_files, readRDS)

#Perturb-seq
pert.seur <- readRDS(paste0(indir.pert, "Seurat_merged_QC_2025-08-04.rds"))
guide_assignment <- readMM(paste0(root, "promoter_screen/interpretation/guide_assignments/matrix/grna_assignment_matrix.mtx"))
guide_rows <- fread(paste0(root, "promoter_screen/interpretation/guide_assignments/matrix/grna_assignment_matrix_rownames.txt"), header = F)
guide_cols <- fread(paste0(root, "promoter_screen/interpretation/guide_assignments/matrix/grna_assignment_matrix_colnames.txt"), header = F)

#enhancer screen results 
enh <- fread(paste0(enrdir, "input/cre_perturbation_screen_results_traits.csv"), )

#genes we tested
genes_to_test <- fread("/g/stegle/schrod/code/TCell/gene_set.csv") #set of genes to test based on being highly variable + our targets


#### ADJUST ####

#keep only relevant enhancer info
enh <- enh %>% 
  filter(is_target == TRUE) %>%
  select(grna_target, response_id, element_type, interaction_type, enhancer_like_interaction) %>%
  distinct()

#### NORMALIZE COUNTS (CP10K + log1p) ####

normalize_log1p <- function(mat, scale_factor = 1e4) {
  # mat: sparse dgCMatrix (genes x cells)
  lib_size <- Matrix::colSums(mat)
  lib_size[lib_size == 0] <- 1  # avoid division by zero
  
  mat_norm <- t( t(mat) / lib_size ) * scale_factor
  mat_norm <- log1p(mat_norm)
  
  return(mat_norm)
}

#Extract count matrices while keeping only the genes that we tested
#Extract gRNA info

#### TAP-seq #### 

## ---- extract + subset gene counts ----
res_list <- lapply(seq_along(tap_list), function(i) {
  
  message("working on panel ", i)
  
  ##counts
  mat <- GetAssayData(tap_list[[i]], assay = "RNA", layer = "counts")
  # keep only tested genes present in this object
  genes_present <- intersect(genes_to_test$bpcells_name, rownames(mat))
  mat <- mat[genes_present, , drop = FALSE]
  # make cell IDs unique
  new_cells <- paste0("panel", i, "_", colnames(mat))
  colnames(mat) <- new_cells
  
  ##gRNA assignments
  meta_dt <- tap_list[[i]]@meta.data %>%
    rownames_to_column(var = "cell_name") %>%
    select(cell = cell_name, gRNAs) %>%
    as.data.table()
  
  meta_dt[, cell := paste0("panel", i, "_", cell)]
  meta_dt <- meta_dt[cell %in% new_cells, .(cell, gRNAs)]
  
  list(counts = mat, grna = meta_dt)
})

## ---- split back out ----
count_list <- lapply(res_list, `[[`, "counts")
grna_list  <- lapply(res_list, `[[`, "grna")

## ---- union only across tested genes ----
all_genes <- Reduce(union, lapply(count_list, rownames))

## ---- pad missing tested genes only ----
pad_counts <- function(mat, all_genes) {
  missing <- setdiff(all_genes, rownames(mat))
  if (length(missing) > 0) {
    zero_block <- Matrix(
      0,
      nrow = length(missing),
      ncol = ncol(mat),
      sparse = TRUE,
      dimnames = list(missing, colnames(mat))
    )
    mat <- rbind(mat, zero_block)
  }
  mat[all_genes, , drop = FALSE]
}

count_list_padded <- lapply(count_list, pad_counts, all_genes = all_genes)
tap_counts <- Reduce(Matrix::cbind2, count_list_padded)

#normalize with log1p
tap_counts  <- normalize_log1p(tap_counts)

#Extract gRNA info
tap_grna <- rbindlist(grna_list)

#Add enhancer annotation
tap_grna <- tap_grna %>%
  mutate(grna_target = sub("-[^-]+$", "", gRNAs)) %>%
  left_join(enh, by = "grna_target", relationship = "many-to-many")

#### Perturb-seq ####
pert_counts <- pert.seur@assays$RNA$counts[genes_to_test$bpcells_name, ]
#normalize with log1p
pert_counts <- normalize_log1p(pert_counts)

# Perturb-seq 
ga_summary <- summary(guide_assignment)
pert_grna <- data.table(
  cell = guide_cols$V1[ga_summary$j],
  guide = guide_rows$V1[ga_summary$i])

#Remove heavy seurat objects
rm(res_list)
rm(tap_list)
rm(pert.seur)
gc()

#### FUNCTION TO EXTRACT DATA FOR A SINGLE GENE ####
get_gene_dt <- function(gene_name, tap_counts, pert_counts, tap_grna, pert_grna) {
  
  ## TAP-seq ##
  tap_counts_gene <- tap_counts[gene_name, , drop = FALSE]  # single gene
  tap_counts_dt <- as.data.table(t(tap_counts_gene), keep.rownames = "cell")
  
  tap_dt <- merge(tap_grna, tap_counts_dt, by = "cell", all = FALSE) %>%
    rename(exp_count = !!gene_name) %>%
    group_by(gRNAs) %>%
    mutate(any_on_target = any(str_detect(response_id, gene_name)),
           guide_type = case_when(any_on_target ~ "on target",
                                  str_detect(gRNAs, "Non-targeting") ~ "non targeting",
                                  TRUE ~ "other"),
           dataset = "TAP-seq") %>%
    select(cell, gRNAs, exp_count, guide_type, dataset, response_id, interaction_type, enhancer_like_interaction)
  
  ## Perturb-seq ##
  pert_counts_gene <- pert_counts[gene_name, , drop = FALSE]
  pert_counts_dt <- as.data.table(
    as.matrix(t(pert_counts_gene)),
    keep.rownames = "cell")
  
  pert_dt <- merge(pert_grna, pert_counts_dt, by = "cell", all = FALSE) %>%
    rename(exp_count = !!gene_name) %>%
    mutate(guide_type = case_when(str_detect(guide, gene_name) ~ "on target",
                                  str_detect(guide, "NO-TARGET") ~ "non targeting",
                                  TRUE ~ "other"),
           dataset = "Perturb-seq",
           response_id = str_extract(guide, "(?<=GUIDE-).*-(?=\\d+$)"),
           response_id = str_remove(response_id, "-$"),  
           interaction_type = "perturb-seq",
           enhancer_like_interaction = FALSE)
  
  ## Combine ##
  combined_dt <- bind_rows(
    tap_dt %>% select(cell, guide = gRNAs, exp_count, guide_type, dataset, response_id, interaction_type, enhancer_like_interaction),
    pert_dt %>% select(cell, guide, exp_count, guide_type, dataset, response_id, interaction_type, enhancer_like_interaction)) %>%
    filter(guide_type != "other")
  
  return(combined_dt)
}

#### PLOT FOR ONE GENE ####

## colors
class_colors <- c(
  "Perturb-seq NT" = "grey50",  
  "Perturb-seq on target" = "goldenrod",  
  "TAP-seq NT" = "grey80",  
  "TAP-seq promoter-proximal" = "#0263a3",  
  "TAP-seq promoter-distal" = "#4299d4", 
  "TAP-seq intragenic-proximal" ="#9ad8fc",
  "TAP-seq intragenic-distal" = "#e57373", 
  "TAP-seq intergenic-distal"= "#c62828" 
)

#select gene
gene <-  "IL2RA" # "IL2RA" | "LTB" | "AAMP" | "NRROS" | "PHF19" | "RBM17"

#get gene expression levels for this gene from perturb- and tap-seq datasets
gene_dt <- get_gene_dt(gene, tap_counts, pert_counts, tap_grna, pert_grna)

#filter and adjust
gene_dt_filt <- gene_dt %>%
  filter(response_id == gene | is.na(response_id) | response_id == "NO-TARGET") %>% #keep only cells with guides whose target was the current gene
  mutate(grna_target = sub("-[^-]+$", "", guide), #save peak name
         panel = str_extract(cell, "^panel\\d+"), #extract panel from cell name for TAP-seq data
         combined_type = case_when(dataset == "Perturb-seq" & guide_type == "non targeting" ~ "Perturb-seq NT",
                                   dataset == "Perturb-seq" & guide_type == "on target" ~ "Perturb-seq on target",
                                   dataset == "TAP-seq" & guide_type == "non targeting" ~ "TAP-seq NT",
                                   dataset == "TAP-seq" & interaction_type == "promoter-proximal" ~ "TAP-seq promoter-proximal",
                                   dataset == "TAP-seq" & interaction_type == "promoter-distal" ~ "TAP-seq promoter-distal",
                                   dataset == "TAP-seq" & interaction_type == "intragenic-proximal" ~ "TAP-seq intragenic-proximal",
                                   dataset == "TAP-seq" & interaction_type == "intragenic-distal" ~ "TAP-seq intragenic-distal",
                                   dataset == "TAP-seq" & interaction_type == "intergenic-distal" ~ "TAP-seq intergenic-distal"),
         plot_peak = case_when(dataset == "Perturb-seq" ~ combined_type,
                          combined_type == "TAP-seq NT" ~ combined_type,
                          TRUE ~ grna_target), #plot per peak
         # plot_guide = case_when(dataset == "Perturb-seq" ~ combined_type,
         #                        combined_type == "TAP-seq NT" ~ combined_type,
         #                        TRUE ~ guide)) #plot per individual guide
         plot_guide = case_when(combined_type == "Perturb-seq NT" ~ combined_type,
                                combined_type == "TAP-seq NT" ~ combined_type, 
                                TRUE ~ guide)) #plot per individual guide

#identify the panel that the current gene was expressed in
#NB DO NOT DO THIS FOR PROMOTER TARGETS READ OUT IN EACH PANEL
panel_nr <- gene_dt_filt %>%
  group_by(panel) %>%
  summarise(n = n()) %>%
  filter(!is.na(panel)) %>%
  slice_max(n) %>%
  pull(panel)

#keep only cells from Perturb-seq or the right panel from TAP-seq to avoid overinflation of 0 values from panels where the gene wasn't read out
gene_dt_filt <- gene_dt_filt %>%
  filter(panel == panel_nr | is.na(panel))

#compute median expression per peak
order_df_peak <- gene_dt_filt %>%
  group_by(combined_type, plot_peak) %>%
  summarise(median_exp = median(exp_count), .groups = "drop") %>%
  arrange(combined_type, median_exp)

#compute median expression per guide
order_df_guide <- gene_dt_filt %>%
  group_by(combined_type, plot_peak, plot_guide) %>%
  summarise(median_exp = median(exp_count), .groups = "drop") %>%
  #arrange(combined_type, median_exp)
  arrange(combined_type, plot_peak, median_exp)

#save median expression values for non-targeting controls
median_pert <- order_df_peak %>%
  filter(combined_type == "Perturb-seq NT") %>%
  pull(median_exp)

median_tap <- order_df_peak %>%
  filter(combined_type == "TAP-seq NT") %>%
  pull(median_exp)

#boxplot of expression values for this gene - split per peak
ggplot(gene_dt_filt, aes(x = fct_relevel(plot_peak, order_df_peak$plot_peak), y = exp_count)) +
  geom_boxplot(aes(fill = combined_type), ) +
  geom_hline(yintercept = median_tap, linetype = "dashed", colour = "grey80") +
  geom_hline(yintercept = median_pert, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = class_colors) +
  ylab(paste0(gene, " expression\n")) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.title = element_blank())
ggsave(paste0(outdir, "Boxplot_TAPseq_Perturbseq_", gene, ".pdf"))

#boxplot of expression values for this gene - split per individual guide
ggplot(gene_dt_filt, aes(x = fct_relevel(plot_guide, order_df_guide$plot_guide), y = exp_count)) +
  geom_boxplot(aes(fill = combined_type), ) +
  geom_hline(yintercept = median_tap, linetype = "dashed", colour = "grey80") +
  geom_hline(yintercept = median_pert, linetype = "dashed", colour = "grey50") +
  scale_fill_manual(values = class_colors) +
  ylab(paste0(gene, " expression\n")) +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 90),
        legend.title = element_blank())
#ggsave(paste0(outdir, "Boxplot_TAPseq_Perturbseq_", gene, "_individual_guides.pdf"), height = 5, width = 14)
#ggsave(paste0(outdir, "Boxplot_TAPseq_Perturbseq_", gene, "_individual_guides_reordered.pdf"), height = 5, width = 12)
ggsave(paste0(outdir, "Boxplot_TAPseq_Perturbseq_", gene, "_individual_guides_CRE_and_promoter.pdf"), height = 5, width = 16)

 




# ggplot(gene_dt_filt, aes(x = guide_type, y = exp_count, fill = dataset)) +
#   geom_hline(yintercept = mean_pert, linetype = "dashed", colour = "grey50") +
#   geom_hline(yintercept = mean_tap, linetype = "dashed", colour = "grey80") +
#   geom_violin(trim = TRUE, fill = "#0263a3") +
#   stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +
#   facet_wrap(~ dataset) +
#   ggtitle(paste0(gene, " expression")) +
#   theme_classic()
# 
# ggplot(gene_dt_filt, aes(x = combined_type, y = exp_count, fill = dataset)) +
#   geom_hline(yintercept = mean_pert, linetype = "dashed", colour = "grey50") +
#   geom_hline(yintercept = mean_tap, linetype = "dashed", colour = "grey80") +
#   geom_violin(trim = TRUE, fill = "#0263a3") +
#   stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +
#   ggtitle(paste0(gene, " expression")) +
#   theme_classic() +
#   theme(axis.title.x = element_blank())
# 
# ggplot(gene_dt_filt, aes(x = combined_type, y = log10(exp_count), fill = dataset)) +
#   #geom_hline(yintercept = log10(mean_pert), linetype = "dashed", colour = "grey50") + #NOT RIGHT HEIGHT
#   #geom_hline(yintercept = log10(mean_tap), linetype = "dashed", colour = "grey80") +
#   geom_violin(trim = TRUE, fill = "#0263a3") +
#   stat_summary(fun = mean, geom = "point", shape = 23, size = 3, fill = "white", color = "black") +
#   ggtitle(paste0(gene, " expression")) +
#   theme_classic() +
#   theme(axis.title.x = element_blank())
# 


# 
# 
# #keep only expression levels for these 6 genes:
# #LTB, AAMP, IL2RA, NRROS, PHF19, RBM17
# #genes_to_keep <- c("LTB", "AAMP", "IL2RA", "NRROS", "PHF19", "RBM17")
# 
# #TAP-SEQ
# #save count matrix
# tap.counts <- tap.seur@assays$RNA$counts
# 
# #filter for genes to keep
# tap.counts.filt <- tap.counts[genes_to_keep[1:3],]
# 
# #change to data table
# tap.counts_dt <- as.data.table(
#   as.matrix(t(tap.counts.filt)),
#   keep.rownames = "cell")
# 
# #save gRNA info
# ga_dt <- tap.seur@meta.data %>%
#   rownames_to_column(var = "cell_name") %>%
#   select(cell = cell_name, gRNAs) %>%
#   as.data.table()
# 
# #merge gRNA info with expression
# tap_dt <- merge(
#   ga_dt,
#   tap.counts_dt,
#   by = "cell",
#   all.y = TRUE)
# 
# #PERTURB-SEQ
# #save count matrix 
# pert.counts <- pert.seur@assays$RNA$counts
# 
# #filter for genes to keep
# pert.counts.filt <- pert.counts[genes_to_keep,]
# 
# #change to data table
# pert.counts_dt <- as.data.table(
#   as.matrix(t(pert.counts.filt)),
#   keep.rownames = "cell")
# 
# #save gRNA info
# ga_summary <- summary(guide_assignment)
# ga_dt <- data.table(
#   guide = guide_rows$V1[ga_summary$i],
#   cell  = guide_cols$V1[ga_summary$j])
# 
# #merge gRNA info with expression
# pert_dt <- merge(
#   ga_dt,
#   pert.counts_dt,
#   by = "cell",
#   all.y = TRUE)
# 
# #remove unused things
# rm(pert.seur)
# rm(pert.counts)
# rm(pert.counts_dt)
# rm(pert.counts.filt)
# rm(guide_assignment)
# 
# 
# #### PLOT ####
# 
# gene <- "IL2RA"
# 
# pert_subset <- pert_dt %>%
#   select(cell, guide, exp_count = !!sym(gene)) %>%
#   filter(str_detect(guide, gene) | str_detect(guide, "NO-TARGET")) %>%
#   mutate(guide_type = case_when(str_detect(guide, gene) ~ "on target",
#                                 str_detect(guide, "NO-TARGET") ~ "non targeting"))
# 
# ggplot(pert_subset) +
#   geom_violin(aes(x = guide_type, y = log10(exp_count))) +
#   ggtitle(paste0(gene, " expression")) +
#   theme_classic()
# 
# tap_subset <- tap_dt %>%
#   select(cell, gRNAs, exp_count = !!sym(gene)) %>%
#   filter(str_detect(gRNAs, "Non-targeting")) %>%
#   mutate(guide_type = case_when(str_detect(gRNAs, gene) ~ "on target",
#                                 str_detect(gRNAs, "Non-targeting") ~ "non targeting"))
# 
# ggplot(tap_subset) +
#   geom_violin(aes(x = guide_type, y = log10(exp_count))) +
#   ggtitle(paste0(gene, " expression")) +
#   theme_classic()
