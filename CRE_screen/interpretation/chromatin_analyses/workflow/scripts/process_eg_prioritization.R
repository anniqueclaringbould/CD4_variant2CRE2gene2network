## Process E-G pair prioritization table for merging with screen results table (based on Annique's
## code)

# save.image("RDA/process_eg_prioritization.rda")
# stop()

suppressPackageStartupMessages({
  library(data.table)
  library(tidyverse)
})

# load information on prioritized E-G pairs
prioritization <- fread(snakemake@input[[1]])

# add prioritization information to chosen eg pairs
prioritization <- prioritization %>%
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

# simplify prioritization sources
prioritization <- prioritization %>% 
  mutate(source_simplified = case_when(str_detect(source, "GRN") ~ "GRN",
                                       str_detect(source, "eQTL") ~ "eQTL",
                                       str_detect(source, "100kb") ~ "100kb window",
                                       str_detect(source, "overlap") ~ "100kb window",
                                       str_detect(source, "stream") ~ "up- and downstream",
                                       str_detect(source, "nonGWAS_Nila_validation") ~ NA, #too few hits to merit their own category
                                       .default = source))

# keep only pair unique identifier and reason for the prioritization, and collapse per pair
prioritization <- prioritization %>%
  select(pair_uid = link_id, prioritization_source = source_simplified) %>%
  distinct() %>% 
  arrange(pair_uid, prioritization_source) %>% 
  group_by(pair_uid) %>% 
  summarize(prioritization_sources = paste(prioritization_source, collapse = ", "))

# save to output file
write_csv(prioritization, file = snakemake@output[[1]])



