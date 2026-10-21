## plot chromatin and regulatory landscape of alternative CRE screen filtering scenarios

# save.image("RDA/plot_filter_scenarios.rda")
# stop()

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
})

# get input files for all filter scenarios, named by filter scenario
pos_files <- unlist(snakemake@input$positives)
neg_files <- unlist(snakemake@input$negatives)
names(pos_files) <- basename(dirname(pos_files))
names(neg_files) <- basename(dirname(neg_files))

# load all filter scenarios including readable names for plots and retain scenarions to plot
all_scenarios <- read_tsv(snakemake@input$scenarios, show_col_types = FALSE) %>%
  filter(filter_scenario %in% names(pos_files))

# load all input files and combine into 2 tables, one for positive hits and one for negatives
positives <- pos_files %>% 
  lapply(FUN = read_tsv, show_col_types = FALSE) %>% 
  bind_rows(.id = "filter_scenario")

negatives <- neg_files %>% 
  lapply(FUN = read_tsv, show_col_types = FALSE) %>% 
  bind_rows(.id = "filter_scenario")

# combine into one table containing positive and negative cres
cres <- bind_rows(positives, negatives)

# features to include in plots
plot_features <- c("log_2_fold_change", "enh_signal_count_DNase-seq", "enh_signal_avg_H3K27ac",
                   "enh_signal_avg_H3K4me1", "enh_signal_avg_H3K4me3")

# reformat for plotting
plot_cres <- cres %>% 
  pivot_longer(cols = all_of(plot_features), names_to = "column_name", values_to = "signal") %>% 
  mutate(column_name = sub("enh_signal_", "", column_name)) %>% 
  mutate(column_name = if_else(column_name == "log_2_fold_change", "LFC", column_name)) %>% 
  separate(column_name, into = c("metric", "assay"), sep = "_", fill = "right") %>% 
  mutate(assay = if_else(metric == "LFC", "Effect size", assay))

# reformat category and make a factor
collapse_cat <- grep("-10k", unique(plot_cres$category), value = TRUE)
cat_levels <- c("hit TSS", "hit < 10k", "hit 10k-100k", "hit 100k-1Mb", "hit >1Mb", "NS")
plot_cres <- plot_cres %>% 
  mutate(category = if_else(category %in% collapse_cat, true = "hit < 10k", false = category)) %>% 
  mutate(category = factor(category, levels = cat_levels))

# reformat filter scenarios and make a factor
plot_cres <- plot_cres %>% 
  left_join(all_scenarios, by = "filter_scenario") %>% 
  mutate(readable_scenario = factor(readable_scenario, levels = all_scenarios$readable_scenario))

## Create plots ------------------------------------------------------------------------------------

# compute sample sizes per category in each filtering scenario (readable_scenario needs to be the
# same factor as in plot_cres, so that labels are dodged in the same order as boxplots)
sample_sizes <- plot_cres %>%
  select(readable_scenario, grna_target, category) %>%
  distinct() %>%
  count(readable_scenario, category, .drop = FALSE)

# colors to use for the different filter scenarios
filter_colors <- deframe(select(all_scenarios, readable_scenario, color))

# make effect size plot
es_plot <- ggplot(filter(plot_cres, assay == "Effect size"),
                  aes(x = category, y = signal, fill = readable_scenario)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, position = position_dodge(width = 0.75)) +
  geom_text(data = sample_sizes,
            aes(x = category, y = 2.6, label = n, group = readable_scenario),
            position = position_dodge(width = 0.75), size = 2.5, vjust = 0, inherit.aes = FALSE,
            angle = 45) +
  labs(title = "Effect size", y = expression("log"[2] ~ "fold change")) +
  scale_fill_manual(values = filter_colors) +
  scale_y_continuous(limits = c(-4, 3)) +
  theme_bw() +
  guides(fill = guide_legend(ncol = 3)) +
  theme(legend.position = "bottom", axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# make DNase-seq signal plot
dnase_plot <- ggplot(filter(plot_cres, assay == "DNase-seq"),
                     aes(x = category, y = signal + 1, fill = readable_scenario)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, position = position_dodge(width = 0.75)) +
  geom_text(data = sample_sizes,
            aes(x = category, y = 20000, label = n, group = readable_scenario),
            position = position_dodge(width = 0.75), size = 2.5, vjust = 0, inherit.aes = FALSE,
            angle = 45) +
  labs(title = "DNase-seq", y = "Read-depth normalized\nsignal") +
  scale_fill_manual(values = filter_colors) +
  scale_y_log10(limits = c(1, 30000)) +
  theme_bw() +
  guides(fill = guide_legend(ncol = 3)) +
  theme(legend.position = "bottom", axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# make H3K27ac signal plot
h3k27ac_plot <- ggplot(filter(plot_cres, assay == "H3K27ac"),
                       aes(x = category, y = signal + 0.01, fill = readable_scenario)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, position = position_dodge(width = 0.75)) +
  geom_text(data = sample_sizes,
            aes(x = category, y = 200, label = n, group = readable_scenario),
            position = position_dodge(width = 0.75), size = 2.5, vjust = 0, inherit.aes = FALSE,
            angle = 45) +
  labs(title = "H3K27ac", y = "FC over control") +
  scale_fill_manual(values = filter_colors) +
  scale_y_log10(limits = c(NA, 400)) +
  theme_bw() +
  guides(fill = guide_legend(ncol = 3)) +
  theme(legend.position = "bottom", axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# make H3K4me1 signal plot
h3k4me1_plot <- ggplot(filter(plot_cres, assay == "H3K4me1"),
                       aes(x = category, y = signal + 0.01, fill = readable_scenario)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, position = position_dodge(width = 0.75)) +
  geom_text(data = sample_sizes,
            aes(x = category, y = 15, label = n, group = readable_scenario),
            position = position_dodge(width = 0.75), size = 2.5, vjust = 0, inherit.aes = FALSE,
            angle = 45) +
  labs(title = "H3K4me1", y = "FC over control") +
  scale_fill_manual(values = filter_colors) +
  scale_y_log10(limits = c(NA, 35)) +
  theme_bw() +
  guides(fill = guide_legend(ncol = 3)) +
  theme(legend.position = "bottom", axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# make H3K4me3 signal plot
h3k4me3_plot <- ggplot(filter(plot_cres, assay == "H3K4me3"),
                       aes(x = category, y = signal + 0.01, fill = readable_scenario)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.5, position = position_dodge(width = 0.75)) +
  geom_text(data = sample_sizes,
            aes(x = category, y = 200, label = n, group = readable_scenario),
            position = position_dodge(width = 0.75), size = 2.5, vjust = 0, inherit.aes = FALSE,
            angle = 45) +
  labs(title = "H3K4me3", y = "FC over control", fill = "Filter\nscenario") +
  scale_fill_manual(values = filter_colors) +
  scale_y_log10(limits = c(NA, 450)) +
  theme_bw() +
  guides(fill = guide_legend(ncol = 3)) +
  theme(legend.position = "bottom", axis.title.x = element_blank(), panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

# arrange plots into one figure
filter_plot <- plot_grid(es_plot + theme(legend.position = "none"),
                         dnase_plot + theme(legend.position = "none"),
                         h3k27ac_plot + theme(legend.position = "none"),
                         h3k4me1_plot + theme(legend.position = "none"),
                         h3k4me3_plot,
                         ncol = 1, rel_heights = c(0.75, 0.75, 0.75, 0.75, 1.1))

# save plot to pdf file
ggsave(filename = snakemake@output[[1]], plot = filter_plot, width = 8, height = 18)
