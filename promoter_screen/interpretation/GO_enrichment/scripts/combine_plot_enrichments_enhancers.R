################################################
## Script to plot pathway enrichment analyses ##
########## Heatmaps based on p-values ##########
################################################

library(tidyverse)
library(data.table)
library(ComplexHeatmap)
library(circlize)

select <- dplyr::select

# -------------------------
# Directories and settings
# -------------------------
root <- "/g/steinmetz/project/otar_2063/"
indir <- paste0(root, "promoter_screen/interpretation/enrichment/output/")
outdir <- paste0(root, "promoter_screen/interpretation/enrichment/output/")
analysis <- "enhancers_hop2" # "enhancers" | "enhancers_hop2"
enrichment_type <- "GO"   # "GO" | "KEGG" | "OpenTargets" | "Hallmark" | "ImmuneSigDB"
date <- "2026-02-18" 

# -------------------------
# Read and combine enrichment files
# -------------------------
files <- list.files(path = paste0(indir, analysis),
                    pattern = paste0(date, ".tsv"),
                    recursive = TRUE)

enr_files <- files[grepl(paste0("/", enrichment_type, "_"), files)]
enr_files <- enr_files[!grepl("all_enhancers", enr_files)] #select all enhancers or only enhancer-like enhancers

enr_long <- rbindlist(lapply(enr_files, function(f){
  
  dt <- fread(file.path(indir, analysis, f))
  
  perturbation <- str_extract(f, "^[^/]+(?=_downstream_genes)")
  
  # if grouped enhancer (contains "__"), keep first enhancer and add "cluster"
  perturbation <- ifelse(
    grepl("__", perturbation, fixed = TRUE),
    paste0(strsplit(perturbation, "__", fixed = TRUE)[[1]][1], "_cluster"),
    perturbation
  )

  dt %>%
    transmute(
      perturbation = perturbation,
      ID = ID,
      padj = p.adjust,
      description = Description,
      pathway = paste0(ID, ":", description)
    )
}))

enr_long <- enr_long %>%
  arrange(padj)

fwrite(enr_long, paste0(outdir, analysis, "/Combined_", enrichment_type, "_enrichments.tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

# -------------------------
# Build p-value matrix
# -------------------------
enr_mat <- enr_long %>%
  select(perturbation, pathway, padj) %>%
  pivot_wider(names_from = perturbation, values_from = padj) %>%
  column_to_rownames("pathway") %>%
  as.matrix()

# missing = not significant
enr_mat[is.na(enr_mat)] <- 1

# significance still based on raw p
sig_mat <- enr_mat < 0.05
sig_counts <- rowSums(sig_mat)

# transform for plotting
plot_mat <- -log10(enr_mat)

# cap infinite values
plot_mat[is.infinite(plot_mat)] <- max(plot_mat[is.finite(plot_mat)], na.rm = TRUE)

# -------------------------
# Identify shared and specific pathways
# -------------------------
shared_pathways   <- names(sig_counts[sig_counts > 1])
specific_pathways <- names(sig_counts[sig_counts == 1])
sig_pathways <- names(sig_counts[sig_counts >= 1])


# -------------------------
# Overall heatmap: cluster rows
# -------------------------

#make matrix with enrichments that are significant in at least one enhancer
plot_shared <- plot_mat[sig_pathways, , drop = FALSE]

# Color function
col_fun <- colorRamp2(
  c(0, -log10(0.2), -log10(0.05), ceiling(max(plot_shared, na.rm = TRUE))),
  c("white", "white", "#FDBB84", "#B30000")
)

h_all <- Heatmap(
  plot_shared,
  name = "-log10(adj.p)",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  border = TRUE,
  column_title = paste0("All ", enrichment_type, " enrichments (p < 0.05)")
)

h_all_long <- Heatmap(
  plot_shared,
  name = "-log10(adj.p)",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  border = TRUE,
  column_title = paste0("All ", enrichment_type, " enrichments (p < 0.05)")
)


# -------------------------
# Shared heatmap: cluster rows
# -------------------------
plot_shared <- plot_mat[shared_pathways, , drop = FALSE]

# Color function
col_fun <- colorRamp2(
  c(0, -log10(0.2), -log10(0.05), ceiling(max(plot_shared, na.rm = TRUE))),
  c("white", "white", "#FDBB84", "#B30000")
)

h_shared <- Heatmap(
  plot_shared,
  name = "-log10(adj.p)",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  border = TRUE,
  column_title = paste0("Shared ", enrichment_type, " enrichments (p < 0.05)")
)

# -------------------------
# Specific heatmap: fully clustered
# -------------------------

# only include pathways significant in exactly one perturbation
specific_pathways <- names(sig_counts[sig_counts == 1])

# subset matrix
plot_specific <- plot_mat[specific_pathways, , drop = FALSE]

# color function
col_fun_specific <- colorRamp2(
  c(0, -log10(0.2), -log10(0.05), ceiling(max(plot_specific, na.rm = TRUE))),
  c("white", "white", "#FDBB84", "#B30000")
)

h_specific <- Heatmap(
  plot_specific,
  name = "-log10(adj.p)",
  col = col_fun_specific,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 6),
  column_names_gp = gpar(fontsize = 8),
  border = TRUE,
  column_title = paste0("Specific ", enrichment_type, " enrichments (p < 0.05)")
)

# -------------------------
# Draw heatmaps
# -------------------------
pdf(paste0(outdir, analysis, "/", enrichment_type, "_heatmap_enrichment_p-values_all_", analysis, ".pdf"), height = 6, width = 4)
h_all
dev.off()

pdf(paste0(outdir, analysis, "/", enrichment_type, "_heatmap_enrichment_p-values_all_", analysis, "_names.pdf"), height = 6, width = 6)
h_all_long
dev.off()

pdf(paste0(outdir, analysis, "/", enrichment_type, "_heatmap_enrichment_p-values_shared_", analysis, ".pdf"), height = 6, width = 6)
h_shared
dev.off()

pdf(paste0(outdir, analysis, "/", enrichment_type, "_heatmap_enrichment_p-values_specific_", analysis, ".pdf"), height = 16, width = 6)
h_specific
dev.off()
 
