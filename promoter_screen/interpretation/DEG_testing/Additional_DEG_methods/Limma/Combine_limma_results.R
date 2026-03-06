##############################################
## Script to combine significant DE results ##
###### limma-trend on individual cells #######
##############################################

#### SETUP ####
#libraries
library(dplyr)
library(tidyverse)
library(limma)
library(edgeR)
library(ggpubr)
library(Seurat)
library(BPCells)
library(Matrix)
library(data.table)
library(ggrastr)

#directories
root <- "/g/steinmetz/project/otar_2063/promoter_screen/"
indir <- paste0(root, "interpretation/DE/output/per_gene/")
outdir <- paste0(root, "interpretation/DE/output/combined/")

#### READ DATA ####
file_name <- "DEG_single_cell_limma_ebayes_split_[0-9]_NT_removed_all_covariates_prefiltered_genes_gene"
#"DEG_single_cell_limma_ebayes_split_[0-9]_random_background_cells_gene"
#"DEG_single_cell_limma_ebayes_split_[0-9]_NT_removed_all_covariates_prefiltered_genes_gene"
#"DEG_single_cell_limma_ebayes_split_[0-9]_NT_removed_all_covariates_0.1pct_overall_0.1pct_per_split_gene"
#"DEG_single_cell_limma_ebayes_split_[0-9]_NT_removed_all_covariates_0.1pct_per_split_gene"
#"DEG_single_cell_limma_ebayes_split_[0-9]_NT_removed_all_covariates_5pct_per_split_gene"
#"DEG_single_cell_limma_ebayes_split_[0-9]_NT_removed_all_covariates_005pct_per_split_gene"
#"DEG_single_cell_limma_ebayes_split_[0-9]_NT_removed_all_covariates_005pct_gene"
#"DEG_single_cell_limma_ebayes_split_[0-9]_NT_removed_all_covariates_gene" 
#"DEG_single_cell_limma_ebayes_split_[0-9]_gene" 
#"DEG_single_cell_limma_ebayes_subset_4e\\+05_5pct_target_gene" 
#"DEG_single_cell_limma_ebayes_subset_4e\\+05_5pct_gene"
files <- list.files(path = indir, pattern = file_name)

date <- "2026-01-21" #"2026-01-19" "2026-01-13" "2026-01-12" #"2026-01-08" #"2025-12-12" #"2025-12-11"
date2 <- "2026-01-22" #"2026-01-09" 

res <- rbindlist(
  lapply(files, function(file){
    #get gene name
    gene <- gsub(file_name, "", file)
    gene <- gsub(paste0("_", date, ".tsv"), "", gene)
    #gene <- gsub(paste0("_", date2, ".tsv"), "", gene) #specific to split data

    #read file
    dm <- fread(paste0(indir, file))
  
    #keep only nominally significant results
    #keep effect of perturbation on the direct target independent of p-value
    dm.filt <- dm %>%
      filter(P.Value < 0.1 | V1 == gene) %>%
      mutate(Perturbation = gene) %>%
      dplyr::select(Perturbation, Target = "V1", everything())
    
  }),
  use.names = TRUE, fill = TRUE
)

res_full <- rbindlist(
  lapply(files, function(file){
    #get gene name
    gene <- gsub(file_name, "", file)
    gene <- gsub(paste0("_", date, ".tsv"), "", gene)
    #gene <- gsub(paste0("_", date2, ".tsv"), "", gene) #specific to split data
    
    #read file
    dm <- fread(paste0(indir, file))
    
    #no filtering
    dm.filt <- dm %>%
      mutate(Perturbation = gene) %>%
      dplyr::select(Perturbation, Target = "V1", everything())
    
  }),
  use.names = TRUE, fill = TRUE
)

#### ADJUST ####

res <- res %>%
  mutate(significant = ifelse(adj.P.Val < 0.05, TRUE, FALSE),
         is_target = case_when(Perturbation == Target ~ TRUE,
                               Perturbation != Target ~ FALSE)) %>%
  distinct()

#filename <- gsub("4e\\\\\\+05", "400000", file_name)
#filename <- gsub("split_", "split", file_name)
filename <- gsub("_\\[0-9\\]", "", file_name)
fwrite(res, paste0(outdir, filename, "_combined_p01_", Sys.Date(), ".tsv.gz"))
fwrite(res_full, paste0(outdir, filename, "_combined_", Sys.Date(), ".tsv.gz"))

res_sig <- res %>%
  filter(significant == TRUE) 

res_targ <- res %>%
  filter(is_target == TRUE) 

res_targ %>%
  group_by(significant) %>%
  tally()

res_sig %>%
  group_by(Target) %>%
  tally() %>%
  arrange(-n)

res_sig %>%
  group_by(Perturbation) %>%
  tally() %>%
  arrange(-n)

p <- res %>%
  filter(P.Value < 0.001) %>%
  ggplot(., aes(x = logFC, y = -log10(adj.P.Val), col = significant)) +
  geom_point(size = 0.5, alpha = 0.5) +
  theme_classic()

p1 <- res %>%
  filter(P.Value < 0.001) %>%
  ggplot(., aes(x = logFC, y = -log10(adj.P.Val), col = is_target)) +
  geom_point(size = 0.5, alpha = 0.5) +
  theme_classic()

p2 <- res %>%
  filter(is_target == TRUE) %>%
  ggplot(., aes(x = logFC, y = -log10(adj.P.Val), colour = significant)) +
  geom_point(size = 0.5, alpha = 0.5) +
  theme_classic()

p3 <- res_sig %>%
  ggplot(., aes(x = logFC, y = -log10(adj.P.Val), col = AveExpr)) +
  geom_point(size = 0.5, alpha = 0.5) +
  theme_classic()

pdf(paste0(outdir, filename, "_Volcano_plots.pdf"))
rasterize(p, dpi = 300)
rasterize(p1, dpi = 300)
rasterize(p2, dpi = 300)
rasterize(p3, dpi = 300)
dev.off()


#Compare two analyses
#res_first <- fread(paste0(outdir, "DEG_single_cell_limma_ebayes_split_gene_combined_p01_2026-01-14.tsv.gz"))
#res_NT.cov <- fread(paste0(outdir, "DEG_single_cell_limma_ebayes_split_NT_removed_all_covariates_gene_combined_p01_2026-01-13.tsv.gz"))
#res_NT.cov_redone <- fread(paste0(outdir, "DEG_single_cell_limma_ebayes_split_NT_removed_all_covariates_gene_combined_p01.tsv.gz"))
#res_NT.cov_redone <- fread(paste0(outdir, "DEG_single_cell_limma_ebayes_split_NT_removed_all_covariates_0.1pct_per_split_gene_combined_p01_2026-01-20.tsv.gz"))
res_NT.8k <- fread(paste0(outdir, "DEG_single_cell_limma_ebayes_split_NT_removed_all_covariates_prefiltered_genes_gene_combined_p01_2026-01-21.tsv.gz"))
res_NT.8kshuffled <- fread(paste0(outdir, "DEG_single_cell_limma_ebayes_split_random_background_cells_gene_combined_p01_2026-01-23.tsv.gz"))

#dm <- res_first %>% inner_join(res_first_redone, by = c("Perturbation", "Target"), suffix = c(".first", ".first_redone"))
#dm <- res_first %>% inner_join(res_NT.cov, by = c("Perturbation", "Target"), suffix = c(".first", ".cov"))
#dm <- res_NT.cov %>% inner_join(res_NT.cov_redone, by = c("Perturbation", "Target"), suffix = c(".cov", ".cov_redone"))
#dm <- res_NT.cov_redone %>% inner_join(res_NT.cov_redone2, by = c("Perturbation", "Target"), suffix = c(".cov5", ".cov005"))
dm <- res_NT.8k %>% inner_join(res_NT.8kshuffled, by = c("Perturbation", "Target"), suffix = c(".8k", ".shuffled"))

p1 <- ggplot(dm) + 
  geom_abline(slope = 1) +
  geom_point(aes(x = -log10(adj.P.Val.8k),
                 y = -log10(adj.P.Val.shuffled)), 
             size = 0.5, alpha = 0.1) +
  theme_classic() +
  xlab("-log10(p-value) Limma results") +
  ylab("-log10(p-value) Limma results gRNAs shuffled")

p2 <- ggplot(dm) + 
  geom_abline(slope = 1) +
  geom_hline(yintercept = 0) +
  geom_vline(xintercept = 0) +
  geom_point(aes(x = logFC.8k,
                 y = logFC.shuffled), 
             size = 0.5, alpha = 0.1) +
  theme_bw() +
  xlab("logFC Limma results") +
  ylab("logFC Limma results gRNAs shuffled")

p1.r <- rasterize(p1, dpi = 300)
p2.r <- rasterize(p2, dpi = 300)

pdf(paste0(outdir, "Comparison_DEG_8k_vs_gRNA_shuffled.pdf"))
p1.r
p2.r
dev.off()

