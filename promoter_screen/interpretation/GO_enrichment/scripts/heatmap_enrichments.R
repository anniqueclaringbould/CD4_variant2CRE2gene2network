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
#indir <- paste0(root, "promoter_screen/interpretation/enrichment/output/")
#outdir <- paste0(root, "promoter_screen/interpretation/enrichment/output/")
indir <- paste0(root, "promoter_screen/interpretation/enrichment/output_gene_names_rematched/")
outdir <- paste0(root, "promoter_screen/interpretation/enrichment/output_gene_names_rematched/")
analysis <-"diseases" #"diseases" | "diseases_hop2"
enrichment_type <- "GO"   # "GO" | "KEGG" | "OpenTargets" | "Hallmark" | "ImmuneSigDB"
date <- "2026-02-24" #"2026-02-17" 

# -------------------------
# Read and combine enrichment files
# -------------------------
files <- list.files(path = paste0(indir, analysis),
                    pattern = paste0(date, ".tsv"),
                    recursive = TRUE)

enr_files <- files[grepl(paste0("/", enrichment_type, "_"), files)]
enr_files <- enr_files[!grepl("Nila|eQTL|GRN", enr_files)]

enr_long <- rbindlist(lapply(enr_files, function(f){
  dt <- fread(file.path(indir, analysis, f))
  perturbation <- str_extract(f, "^[^/]+(?=_downstream_genes)")
  dt %>%
    transmute(
      perturbation = str_replace_all(perturbation, "_", " "),
      ID = ID,
      padj = p.adjust,
      description = Description,
      pathway = paste0(ID, ":", description),
      genes = geneID
    )
}))

enr_long <- enr_long %>%
  arrange(padj) %>%
  group_by(ID) %>%
  mutate(GO_term_sign = sum(padj < 0.05))

fwrite(enr_long, paste0(outdir, analysis, "/Combined_", enrichment_type, "_enrichments", "_", Sys.Date(), ".tsv"), sep = "\t", col.names = TRUE, row.names = FALSE)

# -------------------------
# Build p-value matrix
# -------------------------
enr_mat <- enr_long %>%
  ungroup() %>%
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
# Disease group definition
# -------------------------
disease_to_group <- c(
  # Allergic
  "Allergic disease" = "Allergic",
  "Asthma" = "Allergic",
  "Eczema" = "Allergic",
  "Psoriasis" = "Allergic",
  
  # IBD
  "Celiac disease" = "IBD",
  "Crohns disease" = "IBD",
  "Inflammatory bowel disease" = "IBD",
  "Ulcerative colitis" = "IBD",
  
  # Rheumatic 
  "Ankylosing spondylitis" = "Rheumatic",
  "Gout" = "Rheumatic",
  "Rheumatoid arthritis" = "Rheumatic",
  "Systemic lupus erythematosus" = "Rheumatic",
  "Multiple sclerosis" = "Rheumatic",
  
  # Liver
  "Primary biliary cholangitis" = "Liver",
  "Primary sclerosing cholangitis" = "Liver",
  
  # Diabetes
  "Type 1 diabetes" = "Diabetes"
)

# keep only diseases present
disease_order <- intersect(names(disease_to_group), colnames(enr_mat))

# Preserve the order exactly as defined in disease_to_group
disease_order <- names(disease_to_group)
disease_order <- disease_order[disease_order %in% colnames(plot_mat)]

# reorder matrix
plot_mat <- plot_mat[, disease_order, drop = FALSE]

#Add disease groups for heatmap splits
disease_group <- factor(
  disease_to_group[disease_order],
  levels = c("Allergic", "IBD", "Rheumatic", "Liver", "Diabetes")
)

# -------------------------
# Identify shared and specific pathways
# -------------------------
shared_pathways   <- names(sig_counts[sig_counts > 5])
specific_pathways <- names(sig_counts[sig_counts == 1])

# -------------------------
# Shared heatmap: cluster rows
# -------------------------
plot_shared <- plot_mat[shared_pathways, disease_order, drop = FALSE]

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
  row_names_gp = gpar(fontsize = 6),
  show_row_names = FALSE,
  row_split = 6,                 # <-- number of major clusters
  cluster_row_slices = TRUE,     # keep clustering within clusters
  row_title = NULL,
  column_split = disease_group,
  column_gap = unit( 2, "mm"),
  border = TRUE,
  cluster_column_slices = FALSE,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  column_names_gp = gpar(fontsize = 8),
  column_title = paste0("Shared ", enrichment_type, " enrichments (p < 0.05)")
)

h_shared_long <- Heatmap(
  plot_shared,
  name = "-log10(adj.p)",
  col = col_fun,
  cluster_rows = TRUE,
  row_names_gp = gpar(fontsize = 6),
  show_row_names = TRUE,
  row_split = 6,                 # <-- number of major clusters
  cluster_row_slices = TRUE,     # keep clustering within clusters
  row_title = NULL,
  column_split = disease_group,
  column_gap = unit( 2, "mm"),
  border = TRUE,
  cluster_column_slices = FALSE,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  column_names_gp = gpar(fontsize = 8),
  column_title = paste0("Shared ", enrichment_type, " enrichments (p < 0.05)")
)

# -------------------------
# Specific heatmap: split by dominant disease and
# order within disease by enrichment strength
# -------------------------

# only include pathways significant in exactly one disease
specific_pathways <- names(sig_counts[sig_counts == 1])

# subset matrix
plot_specific <- plot_mat[specific_pathways, disease_order, drop = FALSE]

# Color function
col_fun <- colorRamp2(
  c(0, -log10(0.2), -log10(0.05), ceiling(max(plot_specific, na.rm = TRUE))),
  c("white", "white", "#FDBB84", "#B30000")
)

# build dataframe
specific_df <- as.data.frame(plot_specific)
specific_df$pathway <- rownames(specific_df)

# dominant disease (column with max -log10(p))
specific_df$disease <- apply(specific_df[, disease_order, drop = FALSE], 1, function(x) {
  disease_order[which.max(x)]
})

# max enrichment in dominant disease
specific_df$max_logp <- apply(specific_df[, disease_order, drop = FALSE], 1, function(x) max(x))

# make disease a factor in column order
specific_df$disease <- factor(specific_df$disease, levels = disease_order)

# order by disease and then descending enrichment
specific_df <- specific_df %>%
  arrange(disease, desc(max_logp), pathway)

# row split for Heatmap
row_split_specific <- specific_df$disease

# rebuild matrix
plot_specific <- as.matrix(specific_df[, disease_order, drop = FALSE])
rownames(plot_specific) <- specific_df$pathway

#HEATMAP
h_specific <- Heatmap(
  plot_specific,
  name = "-log10(adj.p)",
  col = col_fun,
  cluster_rows = FALSE,
  row_names_gp = gpar(fontsize = 6),
  show_row_names = FALSE,
  row_split = row_split_specific,
  row_title = NULL,
  column_split = disease_group,
  column_gap = unit( 2, "mm"),
  border = TRUE,
  cluster_column_slices = FALSE,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  column_names_gp = gpar(fontsize = 8),
  column_title = paste0("Disease-specific ", enrichment_type, " enrichments (p < 0.05)")
)

h_specific_long <- Heatmap(
  plot_specific,
  name = "-log10(adj.p)",
  col = col_fun,
  cluster_rows = FALSE,
  row_names_gp = gpar(fontsize = 6),
  show_row_names = TRUE,
  row_split = row_split_specific,
  row_title = NULL,
  column_split = disease_group,
  column_gap = unit( 2, "mm"),
  border = TRUE,
  cluster_column_slices = FALSE,
  show_column_names = TRUE,
  cluster_columns = FALSE,
  column_names_gp = gpar(fontsize = 8),
  column_title = paste0("Disease-specific ", enrichment_type, " enrichments (p < 0.05)")
)

# -------------------------
# Draw heatmaps
# -------------------------
pdf(paste0(outdir, analysis, "/", enrichment_type, "_heatmap_enrichment_p-values_shared_specific_", analysis, "_", Sys.Date(), ".pdf"), height = 6, width = 6)
h_shared
h_specific
dev.off()

pdf(paste0(outdir, analysis, "/", enrichment_type, "_heatmap_enrichment_p-values_shared_", analysis, "_", Sys.Date(), "_names.pdf"), height = 20, width = 8)
h_shared_long
dev.off()

pdf(paste0(outdir, analysis, "/", enrichment_type, "_heatmap_enrichment_p-values_specific_", analysis, "_", Sys.Date(), "_names.pdf"), height = 36, width = 8)
h_specific_long
dev.off()

