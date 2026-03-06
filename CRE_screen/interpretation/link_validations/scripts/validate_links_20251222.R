##############################################################
## Script to quantify gene-target prioritisation strategies ##
##############################################################

#### SETUP ####
#libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(rtracklayer)
library(UpSetR)

#directories
root <- "/g/steinmetz/project/otar_2063/enhancer_screen/"
outdir <- paste0(root, "interpretation/link_validations/output/")

#### READ DATA ####

#read prioritised gene list
links <- fread(paste0(root, "design/gene_selection/prioritised_gene_list_100kb_up_downstream_eQTL_GRN_aim1_E2G_nonGWAS_extended_windows_2024-03-01.tsv"))

#read results, combined and annotated
res <- fread(paste0(root, "interpretation/chromatin_analyses/results/results_df_with_promoterC_annotated.csv"))

#### ADJUST ####

#Adjust some info in the links dataframe
links <- links %>%
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
  filter(!is.na(dummy)) %>%
  mutate(source_simplified = case_when(str_detect(source, "GRN") ~ "GRN",
                                       str_detect(source, "eQTL") ~ "eQTL",
                                       str_detect(source, "100kb") ~ "100kb window",
                                       str_detect(source, "overlap") ~ "100kb window",
                                       str_detect(source, "stream") ~ "up- and downstream",
                                       str_detect(source, "nonGWAS_Nila_validation") ~ NA, #too few hits to merit their own category
                                       .default = source))

#keep only guide id, gene id, and reason for the prioritisation
links_source <- links %>%
  select(link_id, source_simplified) %>%
  distinct() 

#adjust results to rename control promoters and enhancers
res <- res %>%
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
  mutate(link_id = paste0(grna_target, "_", gene_id)) #save link IDs

#split gene sources
res_gene_sources <- res %>%
  separate_rows(prioritization_sources, sep = ",\\s*")

#annotate links_source with results
res_minimal <- res %>%
  select(link_id, significant, pass_qc, ctrl_perturbation, ctrl_type, interaction_type) %>%
  group_by(link_id) %>%
  mutate(significant = any(significant)) %>%
  mutate(pass_qc = any(pass_qc)) %>%
  distinct()

links_source_annotated <- links_source %>%
  left_join(res_minimal, by = "link_id")

table(links_source_annotated$significant, links_source_annotated$source_simplified)
table(links_source_annotated$interaction_type, links_source_annotated$source_simplified, links_source_annotated$significant)

#GENE SOURCES
#all that pass QC
total_source <- res_gene_sources %>% 
  filter(pass_qc == T) %>% 
  select(link_id, prioritization_sources) %>% 
  distinct() %>% 
  group_by(prioritization_sources) %>% 
  tally(name = "total")

#significant
sign_source <- res_gene_sources %>% 
  filter(pass_qc == T) %>% 
  filter(significant == T) %>%
  select(link_id, prioritization_sources) %>% 
  distinct() %>% 
  group_by(prioritization_sources) %>% 
  tally(name = "significant")

#percentages significant of total per categry
source <- total_source %>%
  inner_join(sign_source, by = "prioritization_sources") %>%
  mutate(pct_significant = round((significant/total)*100, digits = 3))

#prioritised vs. not prioritised
prio <- res_gene_sources %>% 
  filter(pass_qc == T) %>% 
  mutate(prioritised = ifelse(prioritization_sources == "not prioritized", FALSE, TRUE)) %>%
  select(link_id, significant, prioritised) %>% 
  distinct()

tab_prio <- table(prio$significant, prio$prioritised)
fisher.test(tab_prio)

#percentage non-prioritised combinations that ends up being significant
pct_nonprio <- tab_prio[2,1]/sum(tab_prio[,1])*100

#percentage prioritesd combinations that ends up being significant
pct_prio <- tab_prio[2,2]/sum(tab_prio[,2])*100

#how much more times are prioritised links significant compared to non-prioritised links?
pct_prio/pct_nonprio

#SIGNIFICANT COMBINATIONS BY GENE SOURCE AND INTERACTION CLASSES
#Save number of significant hits by link source, split by interaction types
summary_int_type <- res_gene_sources %>%
  filter(!is.na(prioritization_sources)) %>% #remove hits without simplified source info
  filter(pass_qc == T) %>% #keep only results that pass QC
  filter(significant == T) %>% #keep only significant results
  filter(ctrl_perturbation == F) %>% #remove control elements
  filter(pert_chr == gene_chr) %>% #remove trans (cross-chromosomal) effects
  group_by(prioritization_sources, interaction_type) %>%
  tally() %>%
  group_by(prioritization_sources) %>%
  mutate(tot_sum = sum(n)) %>%
  ungroup() %>%
  mutate(pct = round(n/tot_sum*100, digits = 3))

#order of sources for this plot
source_order <- c("eQTL", "E2G", "100kb window", "GRN", "up- and downstream", "not prioritized")

#change control enhancer/promoter names to their locations; make link ID; relevel the distance bins
interaction_levels <- c("promoter-proximal", "promoter-distal", "intragenic-proximal", "intragenic-distal", "intergenic-distal")

# colors for each class
class_colors <- c(
  "promoter-proximal" = "#0263a3",   
  "promoter-distal" = "#4299d4", 
  "intragenic-proximal" ="#9ad8fc",
  "intragenic-distal" = "#e57373", 
  "intergenic-distal"= "#c62828" 
)

# Stacked bar plot of numbers
ggplot(summary_int_type) + 
  geom_bar(stat="identity",
           aes(fill=fct_rev(fct_relevel(interaction_type, interaction_levels)), 
               y=n, 
               x=fct_reorder(prioritization_sources, tot_sum))) +
  scale_fill_manual(values = class_colors) +
  theme_classic() +
  coord_flip() +
  xlab("") +
  ylab("\nNumber of significant CTPs") +
  guides(fill = guide_legend(title = "", reverse = T)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10)) +
  geom_text(data = source, 
            aes(label = paste0(pct_significant, "%"), x = prioritization_sources, y = 500),
            position = position_dodge(width = 0.9),
            hjust = -0.1,
            vjust = 0.5,
            size = 6) +
  expand_limits(y = 550)  # to make space for labels
ggsave(paste0(outdir, "Number_of_validated_links_source_classes_type_", Sys.Date(), ".pdf"), height = 4, width = 10)

#order of sources for this plot
source_order <- c("eQTL", "E2G", "100kb window", "GRN", "up- and downstream", "not prioritized")

ggplot(summary_int_type, 
       aes(fill=fct_rev(fct_relevel(interaction_type, interaction_levels)), 
           y=n, 
           x=fct_rev(fct_relevel(prioritization_sources, source_order)))) + 
  geom_bar(position="fill", 
           stat="identity") +
  scale_fill_manual(values = class_colors) +
  theme_classic() +
  coord_flip() +
  xlab("") +
  ylab("\nProportion (significant CTPs)") +
  guides(fill = guide_legend(title = "", reverse = T)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10))
ggsave(paste0(outdir, "Proportion_of_validated_links_source_classes_type_", Sys.Date(), ".pdf"), height = 4, width = 10)


#SIGNIFICANT COMBINATIONS BY ENHANCER STATUS
#Save number of significant hits by link source, split by enhancer status

summary_enh <- res_gene_sources %>%
  filter(!is.na(prioritization_sources)) %>% #remove hits without simplified source info
  filter(significant == T) %>% #keep only significant results
  filter(ctrl_perturbation == F) %>% #remove control elements
  filter(pert_chr == gene_chr) %>% #remove trans (cross-chromosomal) effects
  group_by(prioritization_sources, enhancer_like_interaction) %>%
  tally() %>%
  group_by(prioritization_sources) %>%
  mutate(tot_sum = sum(n)) %>%
  ungroup()

# Stacked bar plot of numbers
ggplot(summary_enh) + 
  geom_bar(stat="identity",
           aes(fill=enhancer_like_interaction, 
               y=n, 
               x=fct_reorder(prioritization_sources, tot_sum))) +
  theme_classic() +
  coord_flip() +
  xlab("") +
  ylab("\nNumber of significant CTPs") +
  guides(fill = guide_legend(title = "Enhancer-like interaction", reverse = T)) +
  scale_fill_manual(values = c("grey70", "#0263a3")) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10))
ggsave(paste0(outdir, "Number_of_validated_links_source_enhancer_type_", Sys.Date(), ".pdf"), height = 4, width = 10)

ggplot(summary_enh) + 
  geom_bar(position="fill",
           stat="identity",
           aes(fill=enhancer_like_interaction, 
               y=n, 
               x=fct_reorder(prioritization_sources, tot_sum))) +
  scale_fill_manual(values = c("grey70", "#0263a3")) +
  theme_classic() +
  coord_flip() +
  xlab("") +
  ylab("\nNumber of significant CTPs") +
  guides(fill = guide_legend(title = "Enhancer-like interaction", reverse = T)) +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 10))
ggsave(paste0(outdir, "Proportion_of_validated_links_source_enhancer_type_", Sys.Date(), ".pdf"), height = 4, width = 10)


# UPSET PLOTS OF PRIORITISED LINKS #

# Upset plots for the simplified categories of links
links_upset_simple <- links %>%
  filter(source != "aim1") %>% #remove the genes that were added for aim1, since they have no guide connected
  select(gene_id, peak, link_id, source_simplified, dummy) %>%
  distinct() %>%
  pivot_wider(names_from = source_simplified, values_from = dummy, values_fill = 0) %>%
  as.data.frame() 

# All links
set_order = c("eQTL", "E2G", "GRN", "up- and downstream", "100kb window")

pdf(paste0(outdir, "Upset_plot_all_prioritised_links_simplified_", Sys.Date(), ".pdf"), height = 6, width = 5)
upset(links_upset_simple, 
      sets = set_order,
      keep.order = T,
      set_size.show = T,
      show.numbers = "yes",
      number.angles = 0,
      point.size = 2, line.size = 1, 
      mainbar.y.label = "Source intersections", sets.x.label = "Peak-gene links per source")
dev.off()


#PERCENTAGES
#Percentages mentioned in manuscript

#eQTL links proximal
summary_int_type %>% 
  filter(prioritization_sources == "eQTL") %>% 
  filter(str_detect(interaction_type, "proximal")) %>% 
  summarise(sum = sum(pct))

#GRN links distal
summary_int_type %>% 
  filter(prioritization_sources == "GRN") %>% 
  filter(str_detect(interaction_type, "distal")) %>% 
  summarise(sum = sum(pct))

#source type % of links significant
table(links_source_annotated$interaction_type, links_source_annotated$significant)
distal <- links_source_annotated %>%
  group_by(interaction_type, source_simplified) %>%
  summarize(nr = dplyr::n(),
            nr_sign = sum(significant, na.rm = T),
            .groups = "drop") %>%
  mutate(pct_sign_from_all_prioritised = round(nr_sign/length(unique(links_source_annotated$link_id))*100, 2),
         pct_sign_from_prioritised = round(nr_sign/nr*100, 2)) %>%
  arrange(-pct_sign_from_prioritised) 

ggplot(distal) +
  geom_point(aes(y = nr, x = nr_sign, colour = source_simplified, shape = interaction_type)) +
  theme_classic() +
  ylab("Number of prioritised links") +
  xlab("Number of significant links")
             