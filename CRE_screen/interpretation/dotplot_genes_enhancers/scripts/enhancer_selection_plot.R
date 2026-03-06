###########################################################
## Script to plot enhancer activities used for selection ##
###########################################################

#### SETUP ####
#libraries
library(data.table)
library(tidyverse)
library(ggplot2)

#directories
root <- "/g/steinmetz/project/otar/"
res_dir <- "/g/steinmetz/moonen/Rstudio/screen_analyses_DM/results/"
outdir <- paste0(root, "interpretation/dotplot_genes_enhancers/output/")

#### READ DATA ####

#read results, combined and annotated
res <- fread("/g/steinmetz/gschwind/otar/manuscript_analyses/results/results_df_with_promoterC_annotated_chromFeatures.csv")

#read enhancer input data (unfiltered)
enhs <- fread(paste0(root, "select_enhancer_peaks/normalised_peaks_across_conditions/ATAC_counts_normToMax_quantileNorm_euclideanNorm.txt"))


#### ADJUST ####

#change grna target from names to coordinates for controls
#keep only unique enhancer names to ensure we colour only those that were actually tested
res_enhancers <- res %>%
  mutate(grna_target = case_when(grna_target == "CD28" ~ "chr2:203706639-203706639", #add promoter controls as gene names (as in results)
                                 grna_target == "CTLA4" ~ "chr2:203867771-203867771",
                                 grna_target == "IL2RA" ~ "chr10:6062367-6062367",
                                 grna_target == "CD81" ~ "chr11:2377310-2377310",
                                 grna_target == "CD4" ~ "chr12:6789528-6789528",
                                 grna_target == "SP140" ~ "chr2:230225130-230226457", #add pilot controls as gene names (as in results)
                                 grna_target == "CAB39" ~ "chr2:230658624-230659752",
                                 grna_target == "PIM1" ~ "chr6:37050015-37051033",
                                 grna_target == "CCND2" ~ "chr12:4153422-4154015",
                                 grna_target == "JUND" ~ "chr19:18291227-18293782",
                                 .default = grna_target)) %>%
  select(grna_target) %>%
  distinct() %>%
  mutate(tested = TRUE)

#make enhancer peak column and select only relevant columns
enhs_selection <- enhs %>%
  select(chr, start, end, ATAC_naive_16H_TH0.txt, ATAC_naive_D5_TH0.txt, ATAC_naive_16H_UNS.txt) %>%
  mutate(peak = paste0(chr, ":", start, "-", end)) 

#calculate mean values
tmp <- enhs_selection %>%
  pivot_longer(cols = starts_with("ATAC"), names_to = "condition")

mean_16h <- tmp %>% filter(condition == "ATAC_naive_16H_TH0.txt") %>% summarise(mean = mean(value))
mean_5d <- tmp %>% filter(condition == "ATAC_naive_D5_TH0.txt") %>% summarise(mean = mean(value))

#calculate enhancer activity compared to unstimulated base line for plotting
enhs_selection <- enhs_selection %>%
  mutate(peak_above_average = case_when((ATAC_naive_16H_TH0.txt > mean_16h & ATAC_naive_D5_TH0.txt > mean_5d) ~ "16h and 5d",
                                        ATAC_naive_16H_TH0.txt > mean_16h ~ "16h",
                                        ATAC_naive_D5_TH0.txt > mean_5d ~ "5d",
                                        (ATAC_naive_16H_TH0.txt <= mean_16h & ATAC_naive_D5_TH0.txt <= mean_5d) ~ "neither"),
         delta_naive_16H = ATAC_naive_16H_TH0.txt - ATAC_naive_16H_UNS.txt,
         delta_naive_5D = ATAC_naive_D5_TH0.txt - ATAC_naive_16H_UNS.txt,
         selected = case_when((peak_above_average != "neither" & (delta_naive_16H > 0 | delta_naive_5D > 0)) ~ TRUE,
                              .default = FALSE)) %>%
  left_join(res_enhancers, by = c("peak" = "grna_target")) %>%
  mutate(tested = case_when(is.na(tested) ~ FALSE,
                            .default = tested),
         selected_tested = case_when((selected == T & tested == T) ~ "enhancer included in screen",
                                     (selected == T & tested == F) ~ "enhancer active in stimulated condition,\nno GWAS overlap",
                                     selected == F ~ "enhancer not active in stimulated condition"))
  
table(enhs_selection$tested, enhs_selection$selected, useNA = 'always') #sanity check that only the enhancers selected in this way were tested
table(enhs_selection$selected_tested, useNA = 'always')


#### PLOT ####

ggplot(enhs_selection) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  geom_hline(yintercept = 0, colour = "darkgrey") +
  geom_vline(xintercept = 0, colour = "darkgrey") +
  geom_point(aes(x = delta_naive_16H, y = delta_naive_5D, colour = selected_tested, alpha = tested), size = 0.3, alpha = 0.5) +
  geom_point(data = enhs_selection[enhs_selection$selected_tested == "enhancer included in screen",], aes(x = delta_naive_16H, y = delta_naive_5D, colour = selected_tested, alpha = tested), size = 0.3, alpha = 0.5) +
  scale_colour_manual(values = c('pink', 'maroon', 'grey')) +
  scale_alpha_manual(values = c(0.1, 1), guide = "none") +
  xlab("\nEnhancer activity after 16h\ncompared to unstimulated") +
  ylab("Enhancer activity after 5 days\ncompared to unstimulated\n") +
  guides(colour = guide_legend(title = NULL)) +
  coord_fixed() +  # Ensures a square plot area
  scale_x_continuous(limits = range(c(enhs_selection$delta_naive_16H, enhs_selection$delta_naive_5D))) +
  scale_y_continuous(limits = range(c(enhs_selection$delta_naive_16H, enhs_selection$delta_naive_5D))) +
  theme_classic() 
ggsave(paste0(outdir, 'enhancer_regions_selected_by_activity_', Sys.Date(), '.pdf'), height = 6, width = 9, dpi = 400)
