## Make relative expression levels per TAP-seq target gene panel

# save.image("RDA/tapseq_panel_summary_plots.rda")
# stop()

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(qs)
})

# Process input data -------------------------------------------------------------------------------

# load screen results
enh_results <- read_csv(snakemake@input$enh_results, show_col_types = FALSE)

# load grna targets table for each panel
grna_targets <- read_tsv(snakemake@input$grna_targets, show_col_types = FALSE)

# load grna_assignments table
guide_assignments <- read_tsv(snakemake@input$guide_assignments, show_col_types = FALSE)

# load gene TPM table
gene_tpm <- read_csv(snakemake@input$gene_tpm, show_col_types = FALSE)

# load TAP-seq target gene panels
gene_panels <- qread(snakemake@input$gene_panels)

# get target genes per panel and label control genes
gene_panels <- gene_panels %>% 
  bind_rows(., .id = "panel") %>% 
  mutate(ctrl_gene = if_else(!grepl("cluster", cluster), true = TRUE, false = FALSE)) %>% 
  select(panel, gene_id, ctrl_gene) %>% 
  distinct()

# if a gene is both a control and discovery gene in a panel, classify it as control
gene_panels <- gene_panels %>% 
  group_by(panel, gene_id) %>% 
  slice_max(ctrl_gene)

# add gene TPM data
gene_panels <- left_join(gene_panels, gene_tpm, by = "gene_id")

# filter grna targets table for guide that were assigned to cells only if specified
if (snakemake@params$filt_assigned_guides == TRUE) {
  missing_grnas <- setdiff(grna_targets$grna_id, guide_assignments$grna_id)
  message("Filtering out ", length(missing_grnas), " gRNAs not assigned to any cell")
  grna_targets <- filter(grna_targets, !grna_id %in% missing_grnas)
}

# filter gene panels table for genes in enhancer screen results
if (snakemake@params$filt_detected_genes == TRUE) {
  missing_genes <- setdiff(gene_panels$gene, enh_results$response_id)
  message("Filtering out ", length(missing_genes), " genes not included in screen analyses")
  gene_panels <- filter(gene_panels, !gene %in% missing_genes)
}

# Plot number of gRNAs and target genes per panel --------------------------------------------------

# count guides per target gene panel
guides_per_panel <- grna_targets %>% 
  filter(grna_target != "non-targeting") %>% 
  group_by(panel) %>% 
  summarize(disc_grnas = sum(type == "discovery"),
            ctrl_grnas = sum(type %in% c("promoter_control", "enhancer_control")))

# transform to long format for plotting and create new panel label for plots
guides_per_panel <- guides_per_panel %>% 
  pivot_longer(cols = -panel, names_to = "type", values_to = "grnas") %>% 
  mutate(panel = as.integer(sub("Panel", "", panel)),
         panel_label = paste("Panel", panel)) %>% 
  mutate(panel_label = fct_reorder(panel_label, panel, .desc = TRUE))

# plot the number of discovery and control gRNAs per panel
p1 <- ggplot(guides_per_panel, aes(x = grnas, y = panel_label, fill = type)) +
  geom_bar(stat = "identity") +
  labs(x = "gRNAs", y = "Library / TAP-seq panel", title = "gRNAs\nper library panel") +
  scale_fill_manual(values = c(ctrl_grnas = "steelblue", disc_grnas = "gray50")) +
  theme_bw() +
  theme(panel.grid = element_blank())

# count the number of target genes per panel
genes_per_panel <- gene_panels %>% 
  group_by(panel) %>% 
  summarize(disc_genes = sum(ctrl_gene == FALSE),
            ctrl_genes = sum(ctrl_gene == TRUE))

# transform to long format for plotting and create new panel label for plots
genes_per_panel <- genes_per_panel %>% 
  pivot_longer(cols = -panel, names_to = "type", values_to = "genes") %>% 
  mutate(panel = as.integer(sub("Panel", "", panel)),
         panel_label = paste("Panel", panel)) %>% 
  mutate(panel_label = fct_reorder(panel_label, panel, .desc = TRUE))

# plot the number of discovery and control genes per panel
p2 <- ggplot(genes_per_panel, aes(x = genes, y = panel_label, fill = type)) +
  geom_bar(stat = "identity") +
  labs(x = "Target genes", y = "Library / TAP-seq panel",
       title = "TAP-seq target genes\nper TAP-seq panel") +
  scale_fill_manual(values = c(ctrl_genes = "steelblue", disc_genes = "gray50")) +
  theme_bw() +
  theme(panel.grid = element_blank())

# Plot target gene TPM distributions per panel -----------------------------------------------------

# add panel label to gene panels table and relabel controls for overall figure legend
gene_panels <- gene_panels %>% 
  mutate(panel = as.integer(sub("Panel", "", panel)),
         panel_label = paste("Panel", panel)) %>% 
  mutate(panel_label = fct_reorder(panel_label, panel, .desc = TRUE)) %>% 
  mutate(ctrl_gene = if_else(ctrl_gene == TRUE, true = "Control gRNA/gene",
                             false = "Discovery gRNA/gene"))

# order for correct plot layering
gene_panels <- arrange(gene_panels, desc(ctrl_gene))

# calculate cumulative TPM per panel
cumulative_tpm <- gene_panels %>% 
  group_by(panel_label) %>% 
  summarize(tpm = sum(tpm)) %>% 
  mutate(ctrl_gene = "Cumulative TPM")

# plot TPM per target gene panel
p3 <- ggplot(gene_panels, aes(x = tpm, y = panel_label, color = ctrl_gene)) +
  geom_jitter(height = 0.25, size = 0.75) +
  geom_boxplot(outlier.shape = NA, fill = NA, color = "black") +
  geom_point(data = cumulative_tpm, shape = "|", size = 3.5) +
  labs(x = "Transcript-per-million", y = "Library / TAP-seq panel",
       title = "Target gene TPM\nper TAP-seq panel") +
  scale_color_manual(values = c("Control gRNA/gene" = "steelblue",
                               "Discovery gRNA/gene" = "gray50",
                               "Cumulative TPM" = "red3")) +
  scale_x_log10() +
  theme_bw() +
  theme(panel.grid = element_blank(), legend.title = element_blank())

# Arrange plots into figure ------------------------------------------------------------------------
fig <- plot_grid(
  p1 + theme(legend.position = "none"),
  p2 + theme(legend.position = "none", axis.title.y = element_blank(),
             axis.text.y = element_blank(), axis.ticks.y = element_blank()),
  p3 + theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
             axis.ticks.y = element_blank()),
  nrow = 1, rel_widths = c(0.68, 0.5, 0.95)
  )

# save plot to pdf
ggsave(fig, file = snakemake@output[[1]], height = 4, width = 9)
