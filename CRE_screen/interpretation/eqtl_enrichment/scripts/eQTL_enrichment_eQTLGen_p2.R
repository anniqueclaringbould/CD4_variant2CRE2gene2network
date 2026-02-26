########################################
## Script to quantify eQTL enrichment ##
############## eQTLGen p2 ##############
########################################

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

#directories
root <- "/g/steinmetz/project/otar_2063/enhancer_screen/"
eqtldir <- paste0(root, "interpretation/eqtl_enrichment/")
outdir <- paste0(eqtldir, "output/")

#### READ DATA ####

#read prioritised gene list
links <- fread(paste0(root, "design/gene_selection/prioritised_gene_list_100kb_up_downstream_eQTL_GRN_aim1_E2G_nonGWAS_extended_windows_2024-03-01.tsv"))

#read screen results, combined and annotated
#res <- fread("/g/steinmetz/gschwind/otar/manuscript_analyses/results/results_df_with_promoterC_annotated.csv")
res <- fread(paste0(root, "interpretation/chromatin_analyses/results/results_df_with_promoterC_annotated.csv"))

#read eQTLs from eQTLGen
cis_eqtlgen <- fread(paste0(eqtldir, "data/independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping_cis-eQTL_export.txt.gz"))


#### ADJUST ####

#Separate the information about the distinct links
links_source <- links %>%
  mutate(link_id = paste0(peak, "_", gene_id)) %>%
  select(link_id, source) %>%
  distinct() 

#add the link sources and mark not prioritsed effects
res2 <- res %>%
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
  distinct()


#### COMPARE ####

#make results into GRanges object 
gr.res <- res2 %>% 
  makeGRangesFromDataFrame(., seqnames.field="pert_chr", start.field="pert_start", end.field="pert_end", keep.extra.columns=T)
seqlevelsStyle(gr.res) <- "UCSC" #consistent gene formatting

#make eQTLs into GRanges object 
gr.eqtl_cis <- cis_eqtlgen %>% 
  makeGRangesFromDataFrame(., seqnames.field="chromosome", start.field="bp", end.field="bp", keep.extra.columns=T)
seqlevelsStyle(gr.eqtl_cis) <- "UCSC"

#overlap based on enhancer location and SNP location
gr.res_overlap <- join_overlap_inner(gr.res, gr.eqtl_cis, maxgap = 0)

#keep only results where the eQTL gene also matches the result gene
res_overlap <- gr.res_overlap %>%
  as.data.frame() %>%
  filter(phenotype == gene_id) 

#consolidate overlapping results to one line per peak-gene combination
res_overlap_minimal <- res_overlap %>%
  select(link_id, significant, variant) %>%
  distinct() %>%
  group_by(link_id, significant) %>%
  summarise(variant = paste(unique(variant), collapse = ","), .groups = "drop")

table(res_overlap_minimal$significant)

#annotate general results with eQTL overlaps
res_eQTL_ann <- res2 %>%
  left_join(res_overlap_minimal, by = c("link_id", "significant"))
fwrite(res_eQTL_ann, paste0(outdir, "Screen_results_full_annotated_eQTL_", Sys.Date(), ".tsv"), sep = "\t", col.names = T, row.names = F)

#keep one line per enh-gene link
res_eQTL_ann_minimal <- res_eQTL_ann %>%
  select(link_id, response_id, is_target, interaction_type, variant, significant) %>%
  group_by(link_id) %>%
  mutate(significant = any(significant)) %>%
  distinct()
fwrite(res_eQTL_ann_minimal, paste0(outdir, "Screen_results_minimal_annotated_eQTL_", Sys.Date(), ".tsv"), sep = "\t", col.names = T, row.names = F)

#check that worked
nrow(res_eQTL_ann_minimal) == length(unique(res2$link_id))

#table of eQTL annotation of significant vs. non-significant effects
cont.table.sign <- table(is.na(res_eQTL_ann_minimal$variant), res_eQTL_ann_minimal$significant)  
  
fisher.test(cont.table.sign)
chisq.test(cont.table.sign)

# Calculate log odds ratio
sign_eqtl <- cont.table.sign[1,2]
nonsign_eqtl <- cont.table.sign[1,1]
sign_noeqtl <- cont.table.sign[2,2]
nonsign_noeqtl <- cont.table.sign[2,1]

# Calculate Odds Ratio
odds_ratio.sign <- (sign_eqtl * nonsign_noeqtl) / (nonsign_eqtl * sign_noeqtl)
log_odds_ratio.sign <- log(odds_ratio.sign)

#table of eQTL annotation of prioritised vs. non-prioritised effects
cont.table.prio <- table(is.na(res_eQTL_ann_minimal$variant), res_eQTL_ann_minimal$is_target)  

fisher.test(cont.table.prio)
chisq.test(cont.table.prio)

# Calculate log odds ratio
sign_eqtl <- cont.table.prio[1,2]
nonsign_eqtl <- cont.table.prio[1,1]
sign_noeqtl <- cont.table.prio[2,2]
nonsign_noeqtl <- cont.table.prio[2,1]

# Calculate Odds Ratio
odds_ratio.prio <- (sign_eqtl * nonsign_noeqtl) / (nonsign_eqtl * sign_noeqtl)
log_odds_ratio.prio <- log(odds_ratio.prio)

#table of eQTL annotation of significant prioritised vs. non-significant prioritised effects
cont.table.prio.sign <- res_eQTL_ann_minimal %>%
  filter(is_target == TRUE) %>%
  mutate(any_variant = case_when(is.na(variant) ~ FALSE,
                                 !is.na(variant) ~ TRUE)) %>%
  group_by(significant, any_variant) %>%
  tally() %>%
  pivot_wider(names_from = any_variant, values_from = n)

fisher.test(cont.table.prio.sign)
chisq.test(cont.table.prio.sign)

# Calculate log odds ratio
sign_eqtl <- cont.table.prio.sign[2,3]
nonsign_eqtl <- cont.table.prio.sign[1,3]
sign_noeqtl <- cont.table.prio.sign[2,2]
nonsign_noeqtl <- cont.table.prio.sign[1,2]

# Calculate Odds Ratio
odds_ratio.prio.sign <- (sign_eqtl * nonsign_noeqtl) / (nonsign_eqtl * sign_noeqtl)
log_odds_ratio.prio.sign <- log(odds_ratio.prio.sign)

#table of interaction types that were validated by eQTL
res_eQTL_ann_minimal %>%
  filter(!is.na(variant)) %>%
  group_by(interaction_type, significant) %>% 
  tally()

#distal effects validated by eQTL
examples <- res_eQTL_ann_minimal %>%
  filter(!is.na(variant)) %>%
  filter(is_target == TRUE) %>%
  filter(significant == TRUE) %>%
  filter(interaction_type == "distal")
fwrite(examples, paste0(outdir, "Distal_interactions_significant_validated_by_eQTL_", Sys.Date(), ".tsv"), sep = "\t", col.names = T, row.names = F)

#write as bedpe
examples_bedpe <- res_eQTL_ann %>%
  filter(!is.na(variant)) %>%
  filter(is_target == TRUE) %>%
  filter(significant == TRUE) %>%
  filter(interaction_type == "distal") %>%
  select(link_id, response_id, gene_id, gene_chr, gene_start, gene_end, variant, significant) %>%
  group_by(link_id) %>%
  mutate(significant = any(significant)) %>%
  distinct() %>%
  left_join(cis_eqtlgen, by = c("variant", "gene_id" = "phenotype")) %>%
  ungroup() %>%
  mutate(snp_chr = paste0("chr", chromosome),
         snp_end = bp + 1) %>%
  select(snp_chr, snp_start = bp, snp_end, gene_chr, gene_start, gene_end, link_id)
fwrite(examples_bedpe, paste0(outdir, "Distal_interactions_significant_validated_by_eQTL_", Sys.Date(), ".bedpe"), sep = "\t", col.names = F, row.names = F)


#### PLOT ####

#LOG ODDS RATIOS
log_odds <- data.frame(c(log_odds_ratio.prio, log_odds_ratio.sign, log_odds_ratio.prio.sign$`TRUE`),
                       test = c("prioritised", "significant", "prioritised and\nsignificant"))
colnames(log_odds) <- c("log(OR)", "test")

my_colours <- RColorBrewer::brewer.pal(4, "YlGnBu")[2:5]

#Barplot of eQTL enrichment
ggplot(log_odds, aes(fill = test, y=`log(OR)`, x=fct_reorder(test, `log(OR)`))) + 
  geom_bar(stat="identity") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50", linewidth = 1, alpha = 0.5) +
  theme_classic() +
  coord_flip() +
  scale_fill_manual(values = my_colours) +
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 16),
        title = element_text(size = 20),
        legend.position = "none") +
  xlab("") +
  ylab("\nlog(odds ratio) of eQTL enrichment") +
  geom_text(data = log_odds,
            aes(label = test, y = 0.05),
            position = position_dodge(width = 0.9),
            hjust = -0.1,
            vjust = 0.5,
            size = 6,
            colour = "white")
ggsave(paste0(outdir, "Log_odds_ratio_eQTL_enrichment_", Sys.Date(), ".pdf"), height = 3.5, width = 7)


#INTERACTION TYPES
#compare interaction types between significant ETPs with eQTL and without eQTL
eqtl_interaction_types <- res_eQTL_ann_minimal %>%
  filter(is_target == TRUE) %>% #keep only prioritised hits
  mutate(any_variant = case_when(is.na(variant) ~ FALSE,
                                 !is.na(variant) ~ TRUE)) %>%
  group_by(any_variant, interaction_type) %>%
  tally() %>%
  mutate(eQTL = case_when(any_variant == TRUE ~ "eQTL",
                          any_variant == FALSE ~ "no eQTL"))

# Colour palette to match with prior colouring
# colors for each class
class_colors <- c(
  "promoter-proximal" = "#0263a3",   
  "promoter-distal" = "#4299d4", 
  "intragenic-proximal" ="#9ad8fc",
  "intragenic-distal" = "#e57373", 
  "intergenic-distal"= "#c62828" 
)
#Order of interaction types
interaction_levels <- c("promoter-proximal", "promoter-distal", "intragenic-proximal", "intragenic-distal", "intergenic-distal")

# Stacked bar plot of proportions
ggplot(eqtl_interaction_types, aes(fill=fct_rev(fct_relevel(interaction_type, interaction_levels)), y=n, x=fct_rev(eQTL))) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = class_colors)  +
  theme_classic() +
  coord_flip() +
  scale_y_continuous(labels = scales::percent_format(scale = 100)) +
  ylab("\n% of significant prioritised enhancer-gene links") +
  xlab("") +
  guides(fill = guide_legend(title = "", reverse = T))  +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        title = element_text(size = 16),
        legend.text = element_text(size = 12))
ggsave(paste0(outdir, "Percentage_interaction_types_prioritised_effects_eQTL_", Sys.Date(), ".pdf"), height = 2, width = 8)


  
  
