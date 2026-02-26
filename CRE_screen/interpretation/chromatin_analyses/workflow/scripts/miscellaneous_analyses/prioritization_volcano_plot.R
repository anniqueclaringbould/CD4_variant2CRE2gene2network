## Make enhancer screen volcano plot with pairs colored by prioritization source

# save.image("RDA/prioritization_volcano_plot.rda")
# stop()

suppressPackageStartupMessages({
  library(tidyverse)
  library(data.table)
})

# Add prioritization sources to enhancer screen results --------------------------------------------

# load enhancer screen results
results <- read_csv(snakemake@input$enh_results, show_col_types = FALSE)

# remove any pairs no passing QC
results <- filter(results, pass_qc == TRUE)

# load information on prioritized E-G pairs
pr_classes <- fread(snakemake@input$pr_classes)

# add prioritization information to chosen eg pairs
pr_classes <- pr_classes %>%
  mutate(SNP_trait = case_when(SNP_trait == "" ~ source, #add source to SNP trait for non GWAS links
                               .default = SNP_trait),
         peak_type = case_when(peak_type == "" ~ "control", #add peak type to peaks that were selected as controls
                               .default = peak_type)) %>%
  filter(SNP_trait != "nonGWAS_control_gsk") %>% #these are in there double
  filter(SNP_trait != "nonGWAS_control_pilot") %>% #these are in there double
  filter(SNP_trait != "nonGWAS_control_promoter") %>% #these are in there double
  filter(SNP_trait != "nonGWAS_eQTL") %>% #these are in there double
  filter(SNP_trait != "nonGWAS_Nila_validation") %>% #these are in there double
  mutate(link_id = paste0(peak, "_", gene_id),
         dummy = 1) %>% #make dummy variable to spread selection categories
  pivot_wider(names_from = source, values_from = dummy) %>%
  mutate(downstream = case_when(downstream == 1 & `100kb window` == 1 ~ NA, #remove annotations as downstream if the link is also in 100kb or overlap
                                downstream == 1 & overlap == 1 ~ NA,
                                .default = downstream),
         upstream = case_when(upstream == 1 & `100kb window` == 1 ~ NA, #remove annotations as upstream if the link is also in 100kb or overlap
                              upstream == 1 & overlap == 1 ~ NA,
                              .default = upstream),
         `100kb window` = case_when(`100kb window` == 1 & overlap == 1 ~ NA, #remove annotations as 100kb window if the link is also overlapping a gene
                                    .default = `100kb window`)) %>%
  pivot_longer(cols = `100kb window`:aim1,
               names_to = "source",
               values_to = "dummy") %>%
  filter(!is.na(dummy))

# keep only guide id, gene id, and reason for the prioritization
pr_classes_source <- pr_classes %>%
  select(link_id, source) %>%
  distinct() 

# add unique identifiers for pairs to results
results <- results %>%
  mutate(grna_target = case_when(grna_target == "CD28" ~ "chr2:203706639-203706639", # add promoter controls as gene names (as in results)
                                 grna_target == "CTLA4" ~ "chr2:203867771-203867771",
                                 grna_target == "IL2RA" ~ "chr10:6062367-6062367",
                                 grna_target == "CD81" ~ "chr11:2377310-2377310",
                                 grna_target == "CD4" ~ "chr12:6789528-6789528",
                                 grna_target == "SP140" ~ "chr2:230225130-230226457", # add pilot controls as gene names (as in results)
                                 grna_target == "CAB39" ~ "chr2:230658624-230659752",
                                 grna_target == "PIM1" ~ "chr6:37050015-37051033",
                                 grna_target == "CCND2" ~ "chr12:4153422-4154015",
                                 grna_target == "JUND" ~ "chr19:18291227-18293782",
                                 .default = grna_target)) %>%
  mutate(link_id = paste0(grna_target, "_", gene_id))

# add the link sources and add simplified source (combine eQTL sources, combine GRN sources, combine
# up- and downstream, combine 100kb and overlap)
results <- results %>%
  left_join(pr_classes_source, by = "link_id", relationship = 'many-to-many') %>%
  mutate(source = case_when(is.na(source) ~ "not prioritized",
                            .default = source)) %>%
  mutate(source_simplified = case_when(str_detect(source, "GRN") ~ "GRN",
                                       str_detect(source, "eQTL") ~ "eQTL",
                                       str_detect(source, "100kb") ~ "100kb window",
                                       str_detect(source, "overlap") ~ "100kb window",
                                       str_detect(source, "stream") ~ "up- and\ndownstream",
                                       str_detect(source, "nonGWAS_Nila_validation") ~ NA, #too few hits to merit their own category
                                       .default = source))

# Make volcano plot  -------------------------------------------------------------------------------

# order results for proper layering of points
results <- results %>%
  mutate(plot_order = if_else(source_simplified == "not prioritized", true = 1, false = 2)) %>% 
  arrange(significant, plot_order) %>% 
  mutate(plot_label = if_else(significant == TRUE,
                              true = paste("Sig.", source_simplified),
                              false = paste("NS", source_simplified)))

# colors for the different categories
cat_colors <- c("Sig. E2G" = "#e41a1c", "Sig. eQTL" = "#377eb8", "Sig. GRN" = "#4daf4a",
                "Sig. 100kb window" = "#984ea3", "Sig. up- and\ndownstream" = "#de6e00",
                "Sig. not prioritized" = "gray42",
                "NS E2G" = "#ffabac", "NS eQTL" = "#9ec6e6", "NS GRN" = "#a1e09f",
                "NS 100kb window" = "#e6b4ed", "NS up- and\ndownstream" = "#faaf64",
                "NS not prioritized" = "gray80")

# maximum absolute effect size to set x-axis limits
max_lfc <- max(abs(results$log_2_fold_change))

# make plot
p <- ggplot(results, aes(x = log_2_fold_change, y = -log10(p_value), color = plot_label)) +
  geom_point(size = 1) +
  geom_vline(xintercept = 0, color = "black", linewidth = 0.25) +
  labs(x = expression("log"[2] ~ "fold change"), y = expression("-log"[10] ~ "(p-value)"),
       color = "ETP category") +
  scale_color_manual(values = cat_colors) +
  scale_x_continuous(limits = c(-max_lfc, max_lfc)) +
  theme_bw() +
  theme(panel.grid = element_blank())

# save plot to file
ggsave(p, filename = snakemake@output[[1]], height = 4, width = 6)
