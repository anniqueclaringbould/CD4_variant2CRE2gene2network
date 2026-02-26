## Make Hi-C interaction plots

# save.image("RDA/hic_interaction_frequency_plot.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(ggpubr)
})

# load ETP table including Hi-C interaction frequencies
etp <- read_csv(snakemake@input[[1]], show_col_types = FALSE)

# make distance to TSS bin an ordered factor
distToTSS_bin_levels <- c("TSS", "1500-10k", "10k-100k", "100k-1M", "1M-5M")
etp <- mutate(etp, distToTSS_bin = fct_relevel(distToTSS_bin, distToTSS_bin_levels))

# extract E-G pairs to plot
etp <- etp %>% 
  filter(interaction_type %in% snakemake@params$int_types) %>% 
  filter(!distToTSS_bin %in% snakemake@params$remove_dist_bins)

# remove pairs without significance and make label
etp <- etp %>% 
  filter(!is.na(significant)) %>% 
  mutate(sig_label = if_else(significant, true = "Sig.", false = "NS")) %>% 
  add_count(sig_label) %>% 
  mutate(sig_label = paste0(sig_label, " (n=", n, ")")) %>% 
  mutate(sig_label = fct_reorder(sig_label, n)) %>% 
  select(-n)

# create new distance to TSS bin levels which include number of pairs
distToTSS_levels <- etp %>% 
  group_by(distToTSS_bin) %>% 
  summarize(sig = sum(significant == TRUE),
            nonsig = sum(significant == FALSE)) %>% 
  mutate(distToTSS_bin = as.character(distToTSS_bin)) %>% 
  mutate(distToTSS_label = paste0(distToTSS_bin, "\nSig.=", sig, "\nNS=", nonsig)) %>% 
  select(distToTSS_label, distToTSS_bin) %>% 
  deframe()

# change distance to TSS labels in ETP table
etp <- mutate(etp, distToTSS_bin = fct_recode(distToTSS_bin, !!!distToTSS_levels))

# set colors for significance labels which include number of pairs
sig_cols <- structure(c("steelblue", "darkgray"), names = levels(etp$sig_label))

# set NA interaction frequency to 0
etp <- replace_na(etp, replace = list(hic_int_freq = 0))

# calculate distance to TSS in bin increments
etp <- etp %>% 
  mutate(hic_bin_dist = abs(dist_to_tss) %/% snakemake@params$hic_res * snakemake@params$hic_res)

## Make plots --------------------------------------------------------------------------------------

# plot interaction frequencies per distance to TSS bin
ggplot(etp, aes(x = distToTSS_bin, y = hic_int_freq + 1, fill = sig_label)) +
  geom_boxplot() +
  labs(y = "Hi-C interaction frequency + 1", x = "Distance to TSS", fill = "Significance") +
  scale_fill_manual(values = sig_cols) +
  scale_y_log10(limits = c(1, 1000)) +
  theme_bw()

# save to pdf file
ggsave(filename = snakemake@output$main, height = 2.5, width = 3 + 0.33 * n_distinct(etp$distToTSS_bin))
  
# make simple 2 group boxplot
p1 <- ggplot(etp, aes(x = sig_label, y = hic_int_freq + 1, fill = sig_label)) +
  geom_boxplot() +
  labs(y = "Hi-C interaction frequency + 1", x = "") +
  scale_fill_manual(values = c("steelblue", "gray66")) +
  scale_y_log10(limits = c(1, 1000)) +
  theme_bw() +
  theme(legend.position = "none")

# plot interaction frequencies as function of distance to TSS
p2 <- ggplot(etp, aes(x = abs(dist_to_tss) / 1000, y = hic_int_freq + 1, color = sig_label)) +
  geom_point(data = filter(etp, significant == FALSE)) +
  geom_point(data = filter(etp, significant == TRUE)) +
  geom_smooth(data = filter(etp, significant == FALSE), se = FALSE, color = "gray40") +
  geom_smooth(data = filter(etp, significant == TRUE), se = FALSE, color = "steelblue") +
  labs(y = "Hi-C interaction frequency + 1", x = "Distance to TSS (kb)", color = "Significance") +
  scale_color_manual(values = c("gray66", "steelblue")) +
  scale_x_log10(limits = c(1, NA)) +
  scale_y_log10(limits = c(1, 1000)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank())

# arrange into one figure and save to pdf file
plot_grid(p1, p2, nrow = 1, rel_widths = c(0.27, 0.73))
ggsave(filename = snakemake@output$supplementary, height = 3, width = 6)

## Pairwise tests for signficant differences -------------------------------------------------------

# perform individual pairwise tests
pw_tests <- compare_means(hic_int_freq ~ significant, data = etp, group.by = "distToTSS_bin",
                          method = "wilcox.test", p.adjust.method = "holm")

# split category labels in to labels and number of enhancers
pw_tests <- pw_tests %>% 
  select(-c(.y., group1, group2)) %>% 
  separate(distToTSS_bin, into = c("distToTSS_bin", "group1", "n_group1", "group2", "n_group2"),
           sep = "\n|=")

# save to text file
write_csv(pw_tests, file = snakemake@output$pairwise_tests)
