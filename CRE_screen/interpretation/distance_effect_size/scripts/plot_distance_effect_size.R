#################################################################
## Script to make plot of ETP effect sizes and distance to TSS ##
#################################################################

#### SETUP ####
#libraries
library(data.table)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(rtracklayer)
library(UpSetR)
library(GenomicRanges)
library(plyranges)
library(openxlsx)
library(scales)

#directories
root <- "/g/steinmetz/project/otar_2063/enhancer_screen/"
eqtldir <- paste0(root, "interpretation/distance_effect_size/")
outdir <- paste0(eqtldir, "output/")

#### READ DATA ####

#read prioritised gene list
links <- fread(paste0(root, "design/gene_selection/prioritised_gene_list_100kb_up_downstream_eQTL_GRN_aim1_E2G_nonGWAS_extended_windows_2024-03-01.tsv"))

#read screen results, combined and annotated
res <- fread("/g/steinmetz/gschwind/otar/manuscript_analyses/results/results_df_with_promoterC_annotated.csv")

#### ADJUST ####

#Separate the information about the distinct links
links_source <- links %>%
  mutate(link_id = paste0(peak, "_", gene_id)) %>%
  select(link_id, source) %>%
  distinct() 

#add the link sources and mark not prioritsed effects
res_minimal <- res %>%
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
  mutate(link_id = paste0(grna_target, "_", gene_id)) %>%
  left_join(links_source, by = "link_id", relationship = 'many-to-many') %>% #this is where additional rows are created
  mutate(source = case_when(source == "overlap" ~ "100kb window", #combine 100kb window and gene overlap into one category
                            source == "nonGWAS_eQTL" ~ "eQTL", #combine eQTLs into one category
                            str_detect(source, "GRN") ~ "GRN", #combine GRNs into one category
                            str_detect(source, "stream") ~ "up- and downstream", #combine up- and downstream into one category
                            is.na(source) ~ "not_prioritised",
                            .default = source),
         source = str_replace_all(source, "nonGWAS_", "")) %>%
  filter(pass_qc == TRUE) %>% #keep only results that pass QC and therefore were evaluated
  select(link_id, response_id, gene_id, grna_target, is_target, dist_to_tss, log_2_fold_change, p_value, significant, ctrl_perturbation, element_type, interaction_type, distToTSS_bin) %>%
  distinct() %>%
  mutate(label_interaction = case_when(interaction_type == "indirect" ~ "promoter-indirect gene",
                                       interaction_type == "target_gene_promoter" ~ "promoter-direct gene",
                                       interaction_type == "target_gene_body" ~ "genic-direct gene",
                                       TRUE ~ interaction_type))

res_minimal_sign <- res_minimal %>%
  filter(significant == TRUE)

#### PLOT ####

#change control enhancer/promoter names to their locations; make link ID; relevel the distance bins
interaction_levels <- c("promoter-direct gene", "promoter-indirect gene", "genic-direct gene", "distal", "trans")

#log10 scale
signed_log_trans <- function(base = 10) {
  trans_new(
    name = paste0("signed_log-", base),
    transform = function(x) sign(x) * log1p(abs(x)) / log(base),
    inverse = function(x) sign(x) * (base^abs(x) - 1),
    domain = c(-Inf, Inf)
  )
}

ggplot(res_minimal_sign) +
  geom_point(aes(x = dist_to_tss, y = log_2_fold_change, colour = fct_relevel(label_interaction, interaction_levels))) +
  scale_x_continuous(
    trans = signed_log_trans(),
    breaks = c(-1e8, -1e6, -1e4, -1500, 0, 1500, 1e4, 1e6, 1e8),
    labels = c("-100M", "-1M", "-10k", "-1.5k", "0", "1.5k", "10k", "1M", "100M")
  ) +
  geom_vline(xintercept = -1500, linetype = "dashed", aes(size = 1, colour = "navy")) +
  geom_vline(xintercept = 1500, linetype = "dashed", aes(size = 1, colour = "navy")) +
  theme_gray() +
  scale_colour_brewer(palette = "YlGnBu", direction = -1)  +
  theme(axis.text.x = element_text(size = 14),
        title = element_text(size = 16),
        legend.title = element_blank()) +
  xlab("\nDistance to TSS") +
  ylab("\nlog2(fold change)\n") 
ggsave(paste0(outdir, "Distance_effect_size_by_interaction_type_", Sys.Date(), ".pdf"), height = 4, width = 13)

#only show the 1Mb window
ggplot(res_minimal_sign) +
  geom_hline(yintercept = 0, colour = "grey") +
  geom_point(aes(x = dist_to_tss, y = log_2_fold_change, colour = fct_relevel(label_interaction, interaction_levels))) +
  scale_x_continuous(
    trans = signed_log_trans(),
    breaks = c(-1e6, -1e4, -1500, 0, 1500, 1e4, 1e6),
    labels = c("-1M", "-10k", "-1.5k", "0", "1.5k", "10k", "1M"),
    limits=c(-1000000, 1000000)
  ) +
  geom_vline(xintercept = -1500, linetype = "dashed", aes(size = 1, colour = "grey")) +
  geom_vline(xintercept = 1500, linetype = "dashed", aes(size = 1, colour = "grey")) +
  theme_classic() +
  scale_colour_brewer(palette = "YlGnBu", direction = -1)  +
  theme(axis.text.x = element_text(size = 14),
        title = element_text(size = 16),
        legend.title = element_blank()) +
  xlab("\nDistance to TSS") +
  ylab("\nlog2(fold change)\n")
ggsave(paste0(outdir, "Distance_effect_size_by_interaction_type_1Mb_", Sys.Date(), ".pdf"), height = 4, width = 10)
