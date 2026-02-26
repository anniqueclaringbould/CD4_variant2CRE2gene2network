## Make chromatin activity plots

# save.image("RDA/chromatin_assay_plot.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggpubr)
})


## Process input data plots ------------------------------------------------------------------------

# load results table with chromatin signal annotations
results <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

# filter for specified E-G pairs to include in analysis
enh_results <- results %>% 
  filter(!is.na(significant)) %>% 
  filter(interaction_type %in% snakemake@params$int_types) %>% 
  filter(!distToTSS_bin %in% snakemake@params$remove_dist_bins)

# chromatin assays to include in plot
chrom_assays <- c("enh_signal_count_DNase-seq", "enh_signal_avg_H3K27ac", "enh_signal_avg_H3K4me1",
                  "enh_signal_avg_H3K4me3", "enh_signal_CPM_NETCAGE")

# select any elements that were detected as positives and categorize by distance
all_bins <- c("hit TSS", "hit 1500-10k", "hit 10k-100k", "hit 100k-1M", "hit 1M-5M", "NS")
positives <- enh_results %>% 
  filter(significant == TRUE) %>% 
  mutate(category = paste("hit", distToTSS_bin)) %>%
  mutate(category = factor(category, levels = all_bins, ordered = TRUE) )

# select unique enhancers and assays to plot
positives <- positives %>% 
  select(grna_target, significant, category, all_of(chrom_assays)) %>% 
  distinct()

# if an enhancer is counted in multiple categories, assign it to the closest distance to TSS bin
positives <- positives %>% 
  group_by(grna_target) %>% 
  slice_min(category, n = 1)

# filter for enhancers that are not linked to any gene and select columns for plot
negatives <- enh_results %>% 
  filter(significant == FALSE, !grna_target %in% positives$grna_target) %>% 
  mutate(category = factor("NS", levels = all_bins, ordered = TRUE)) %>% 
  select(grna_target, significant, category, all_of(chrom_assays)) %>% 
  distinct()

# combine back into one table and convert to long format for plots
plot_results <- bind_rows(positives, negatives) %>% 
  pivot_longer(cols = all_of(chrom_assays),
               names_to = "column_name", values_to = "signal") %>% 
  mutate(column_name = sub("enh_signal_", "", column_name)) %>% 
  separate(column_name, into = c("metric", "assay"), sep = "_")


## Make plots --------------------------------------------------------------------------------------

# filter out NAs for eRNA signal if specified
if (snakemake@params$remove_erna_zeros == TRUE) {
  plot_results <- filter(plot_results, assay != "NETCAGE" | signal > 0)
}

# perform statistical tests (kruskal-wallis) between categories
kw_tests <- plot_results %>% 
  ungroup() %>% 
  compare_means(formula = signal ~ category, method = "kruskal.test", data = ., group.by = "assay") %>% 
  mutate(p.adj.plot = if_else(p.adj < 0.001, true = "< 0.001",
                              false = paste("=", round(p.adj, digits = 2))))

# get sample sizes per category group
sample_sizes <- plot_results %>% 
  ungroup() %>% 
  select(grna_target, category) %>% 
  distinct() %>% 
  arrange(category) %>% 
  count(category, .drop = FALSE)

# set new labels (factor levels) for each category, which include sample size per category
levels <- paste0(sample_sizes$category, " (n=", sample_sizes$n, ")")
levels(plot_results$category) <- levels

# add mock group to KW test table
kw_tests <- mutate(kw_tests, group1 = levels[[2]], group2 = levels[[3]])

# colors for enhancer categories in figures
category_colors <- structure(c("#daa520", "#0263a3", "#4299d4", "#6ebff5", "#9ad8fc", "#a9a9a9"),
                             names = levels)

# plot DNase-seq accessibility
p1 <- ggplot(filter(plot_results, assay == "DNase-seq"),
             aes(x = category, y = signal + 1, fill = category)) +
  facet_wrap(~ assay, nrow = 1) +
  geom_boxplot() +
  stat_pvalue_manual(data = filter(kw_tests, assay == "DNase-seq"), y.position = 15000, size = 3,
                     label = "P adj. {p.adj.plot}", remove.bracket = TRUE, inherit.aes = FALSE) +
  labs(y = paste0("Read-depth normalized", "\n", "signal")) +
  scale_fill_manual(values = category_colors) +
  scale_y_log10(limits = c(1, 20000)) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# plot ChIP-seq signals
p2 <- ggplot(filter(plot_results, assay %in% c("H3K27ac", "H3K4me1", "H3K4me3")),
             aes(x = category, y = signal + 0.01, fill = category)) +
  facet_wrap(~ assay, nrow = 1) +
  geom_boxplot() +
  stat_pvalue_manual(data = filter(kw_tests, assay %in% c("H3K27ac", "H3K4me1", "H3K4me3")),
                     y.position = 120, size = 3, label = "P adj. {p.adj.plot}",
                     remove.bracket = TRUE, inherit.aes = FALSE) +
  labs(y = "FC over control") +
  scale_fill_manual(values = category_colors) +
  scale_y_log10(limits = c( 0.01, 150)) +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# get eRNA signal data
plot_results_erna <- filter(plot_results, assay == "NETCAGE")
levels(plot_results_erna$category) <- all_bins

# re-calculate sample sizes for eRNA plot
sample_sizes <- plot_results_erna %>% 
  ungroup() %>% 
  select(grna_target, category) %>% 
  distinct() %>% 
  arrange(category) %>% 
  count(category, .drop = FALSE)

# set new labels (factor levels) for each category, which include sample size per category
levels <- paste0(sample_sizes$category, " (n=", sample_sizes$n, ")")
levels(plot_results_erna$category) <- levels

# add mock group to KW test table
kw_tests <- mutate(kw_tests, group1 = levels[[2]], group2 = levels[[3]])

# colors for enhancer categories in figures
category_colors <- structure(c("#daa520", "#0263a3", "#4299d4", "#6ebff5", "#9ad8fc", "#a9a9a9"),
                             names = levels)

# plot eRNA signal
p3 <- ggplot(plot_results_erna,
             aes(x = category, y = log2(signal + 0.001), fill = category)) +
  facet_wrap(~ assay, nrow = 1) +
  geom_boxplot() +
  stat_pvalue_manual(data = filter(kw_tests, assay == "NETCAGE"), y.position = 10, size = 3,
                     label = "P adj. {p.adj.plot}", remove.bracket = TRUE, inherit.aes = FALSE) +
  labs(y = expression("log"[2] ~ "(Counts per million (CPM))")) +
  scale_fill_manual(values = category_colors) +
 # scale_y_log10() +
  theme_bw() +
  theme(legend.position = "none", axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# assemble figure panel and save to pdf
plot <- plot_grid(p1, p2, p3, nrow = 1, rel_widths = c(0.42, 1, 0.42))
ggsave(plot, filename = snakemake@output$plot, height = 3, width = 9)

## Pairwise tests for significant differences -------------------------------------------------------

# perform individual pairwise tests
pw_tests <- compare_means(signal ~ category, group.by = "assay", data = ungroup(plot_results),
                          method = "wilcox.test", p.adjust.method = "holm")

# split category labels in to labels and number of enhancers
pw_tests <- pw_tests %>% 
  separate(group1, into = c("group1", "n_group1"), sep = " \\(n=") %>% 
  separate(group2, into = c("group2", "n_group2"), sep = " \\(n=") %>% 
  mutate(n_group1 = sub("\\)", "", n_group1), n_group2 = sub("\\)", "", n_group2)) %>% 
  select(-.y.)

# save to text file
write_csv(pw_tests, file = snakemake@output$pairwise_tests)
