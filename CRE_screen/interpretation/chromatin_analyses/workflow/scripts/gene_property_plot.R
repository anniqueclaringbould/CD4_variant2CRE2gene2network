
# save.image("RDA/gene_property_plot.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggpubr)
})

# load results table with chromatin signal annotations
results <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

# filter for specified E-G pairs to include in analysis
results <- results %>% 
  filter(!is.na(significant)) %>% 
  filter(interaction_type %in% snakemake@params$int_types) %>% 
  filter(!distToTSS_bin %in% snakemake@params$remove_dist_bins)

# count the number of enhancers per gene
enh_per_genes <- results %>% 
  group_by(gene_id) %>% 
  summarize(positives = sum(significant == TRUE),
            total = n())

# gene properties to include in plot
gene_properties <- c("gene_uniq_expr", "gene_tpm")

# select gene-level data to analyze and add the number of detected enhancers
gene_results <- results %>% 
  select(gene_id, all_of(gene_properties)) %>% 
  distinct() %>% 
  left_join(enh_per_genes, by = "gene_id")

# Make general properties plots --------------------------------------------------------------------

# plot the number of genes per ubiquitous expression class
ggplot(gene_results, aes(x = gene_uniq_expr, fill = gene_uniq_expr)) +
  geom_bar() +
  labs(title = "Genes per expression class", x = "Ubiquitously expressed") +
  theme_bw() +
  theme(legend.position = "none")

ggsave(filename = "WIP/gene_level_plots/genes_per_class.pdf", width = 4, height = 4)

# remove any genes with no ubiquitous expression information
gene_results <- filter(gene_results, !is.na(gene_uniq_expr))

# plot the number of enhancers per gene
ggplot(gene_results, aes(x = positives, fill = gene_uniq_expr)) +
  geom_histogram(binwidth = 1, position = "dodge") +
  labs(x = "Number of enhancers", title = "Enhancers per gene") +
  theme_bw()

# save enhancers-per-gene histogram to file
ggsave(filename = "WIP/gene_level_plots/enhancers_per_gene.pdf", width = 4.5, height = 3)
#ggsave(filename = snakemake@output$histgram, height = 4, width = 4)

# plot the number of enhancers per gene
gene_results %>% 
  filter(positives > 0) %>% 
  ggplot(., aes(x = positives, fill = gene_uniq_expr)) +
    geom_histogram(binwidth = 1, position = "dodge") +
    labs(x = "Number of enhancers", title = "Enhancers per gene (>= 1 hit)") +
    theme_bw()

ggsave(filename = "WIP/gene_level_plots/enhancers_per_gene.pdf", width = 4.5, height = 3)

# Percent positives and expression -----------------------------------------------------------------

# add bins for number of positives
gene_results <- gene_results %>% 
  mutate(positives_bin = case_when(
    positives == 0 ~ "0",
    positives %in% 1:2 ~ "1-2",
    TRUE ~ "3+"
  )) %>% 
  filter(!is.na(gene_uniq_expr))

# count the number of genes per number positive bin and ubiquitous expression class
n_genes_per_class <- gene_results %>%
  group_by(gene_uniq_expr, positives_bin) %>% 
  summarize(genes = n())

# compute relative proportions
n_genes_per_class <- n_genes_per_class %>% 
  mutate(pct_genes = genes / sum(genes))

p1 <- ggplot(n_genes_per_class, aes(x = gene_uniq_expr, y = pct_genes, fill = positives_bin)) +
  geom_bar(stat = "identity") +
  labs(title = "Percent positives", x = "Ubiquitously expressed", y = "% of genes", ) +
  theme_bw() +
  theme(legend.position = "none")

p2 <- ggplot(gene_results, aes(x = gene_uniq_expr, y = gene_tpm, fill = positives_bin)) +
  geom_boxplot() +
  labs(title = "Gene expression", x = "Ubiquitously expressed", y = "Gene TPM bin") +
  scale_y_log10() +
  theme_bw()

plot_grid(p1, p2, rel_widths = c(0.6, 1), nrow = 1)

ggsave(filename = "WIP/gene_level_plots/positives_ubiq_expr.pdf", width = 6.5, height = 3)
#ggsave(filename = snakemake@output$properties, width = 6.5, height = 3)


# Number of enhancers as function of expression ----------------------------------------------------

# plot TPM histograms between ubiquitous and non-ubiquitous genes
ggplot(gene_results, aes(x = gene_tpm, color = gene_uniq_expr)) +
  geom_density() +
#  geom_vline(xintercept = 10) +
  labs(title = "Gene expression distribution", x = "Gene TPM", color = "Ubiq. expr.") +
  scale_x_log10() +
  theme_bw()

# save enhancers-per-gene histogram to file
ggsave(filename = "WIP/gene_level_plots/tpm_distributions.pdf", width = 4, height = 3)

# filter out genes with less than 5 TPM
#gene_results <- filter(gene_results, gene_tpm >= 10)

# add TPM bins
gene_results <- gene_results %>% 
  mutate(gene_tpm_bin = cut(gene_tpm, breaks = c(0, 10, 50, 100, 3500),
                          include.lowest = TRUE))

# calculate percentage of positives per expression and number of positives bin
pos_per_expr_level <- gene_results %>% 
  group_by(gene_uniq_expr, gene_tpm_bin, positives_bin) %>% 
  summarize(genes = n()) %>% 
  mutate(pct_genes = genes / sum(genes))

# plot percentage of genes with at least 1 enhancer
ggplot(filter(pos_per_expr_level, positives_bin != "0"),
       aes(x = gene_tpm_bin, pct_genes, fill = gene_uniq_expr)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Genes with at least 1 distal enhancer", x = "Gene TPM bin", y = "% of genes",
       fill = "Ubiq. expr.") +
  theme_bw()

ggsave(filename = "WIP/gene_level_plots/pct_positives_vs_tpm.pdf", width = 5, height = 3)

# break up by number of enhancers
ggplot(filter(pos_per_expr_level, positives_bin != "0"),
       aes(x = gene_tpm_bin, pct_genes, fill = gene_uniq_expr)) +
  facet_wrap(~positives_bin) + 
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of positive distal enhancer", x = "Gene TPM bin", y = "% of genes",
       fill = "Ubiq. expr.") +
  theme_bw()

ggsave(filename = "WIP/gene_level_plots/pct_positives_vs_tpm_nEnh.pdf", width = 8, height = 3)



# Other properties ---------------------------------------------------------------------------------

# chromatin assays to include in plot
gene_properties <- c("gene_signal_count_DNase-seq", "gene_signal_avg_H3K4me3", "gene_tpm",
                     "gene_constraints")

# select gene-level data to plot
gene_results <- results %>% 
  select(gene_id, all_of(gene_properties)) %>% 
  distinct()

# convert to long format and add number of detected enhancers
gene_results <- gene_results %>% 
  pivot_longer(cols = -gene_id, names_to = "property", values_to = "value") %>% 
  left_join(enh_per_genes, by = "gene_id")

# reformat measured property
gene_results <- gene_results %>% 
  mutate(property = case_when(
    property == "gene_signal_count_DNase-seq" ~ "DNase-seq",
    property == "gene_signal_avg_H3K4me3" ~ "H3K4me3",
    property == "gene_tpm" ~ "Expression",
    property == "gene_constraints" ~ "Perturbability",
    TRUE ~ sub("gene_", "", property)
  ))

# create factor for number of hits
gene_results <- gene_results %>% 
  mutate(positive_bin = case_when(
    positives == 0 ~ "0",
    positives == 1 ~ "1",
    positives < 4 ~ "2-3",
    TRUE ~ "4+"
  )) %>% 
  mutate(positive_bin = factor(positive_bin, levels = c("4+", "2-3", "1", "0"), ordered = TRUE))

# colors for categories in figures
category_colors <- structure(c("#0263a3", "#4299d4", "#6ebff5", "#a9a9a9"),
                             names = c("4+", "2-3", "1", "0"))

# plot expression levels of genes
p1 <- ggplot(filter(gene_results, property == "Expression"),
       aes(x = positive_bin, y = value, fill = positive_bin)) +
  facet_wrap(~ property, nrow = 1) +
  geom_boxplot() +
  labs(y = "TPM") +
  scale_fill_manual(values = category_colors) +
  scale_y_log10() +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank())

# plot "Perturbability" of each gene
p2 <- ggplot(filter(gene_results, property == "Perturbability"),
             aes(x = positive_bin, y = value, fill = positive_bin)) +
  facet_wrap(~ property, nrow = 1) +
  geom_boxplot() +
  labs(y = "LoF pLIl") +
  scale_y_log10() +
  scale_fill_manual(values = category_colors) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank())

# plot DNase-seq accessibility of gene promoters
p3 <- ggplot(filter(gene_results, property == "DNase-seq"),
       aes(x = positive_bin, y = value, fill = positive_bin)) +
  facet_wrap(~ property, nrow = 1) +
  geom_boxplot() +
  labs(y = paste0("Read-depth normalized", "\n", "signal")) +
   scale_fill_manual(values = category_colors) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank())

# plot H3K4me3 signal of gene promoters
p4 <- ggplot(filter(gene_results, property == "H3K4me3"),
             aes(x = positive_bin, y = value, fill = positive_bin)) +
  facet_wrap(~ property, nrow = 1) +
  geom_boxplot() +
  labs(y = "FC over control") +
  scale_fill_manual(values = category_colors) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank())

# assemble figure panel and save to pdf
plot <- plot_grid(p1, p2, p3, p4, nrow = 1, rel_widths = c(1.05, 1.08, 1.1, 1))
#ggsave(plot, filename = snakemake@output$properties, height = 3, width = 9)
