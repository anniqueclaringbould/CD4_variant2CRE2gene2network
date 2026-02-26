########################################
## Script to quantify eQTL enrichment ##
############ eQTL catalog ##############
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
indir <- paste0(eqtldir, "data/eQTL_catalogue/")
outdir <- paste0(eqtldir, "output/eQTL_catalogue/")

#### READ DATA ####

#read prioritised gene list
links <- fread(paste0(root, "design/gene_selection/prioritised_gene_list_100kb_up_downstream_eQTL_GRN_aim1_E2G_nonGWAS_extended_windows_2024-03-01.tsv"))

#read screen results, combined and annotated
#res <- fread("/g/steinmetz/gschwind/otar/manuscript_analyses/results/results_df_with_promoterC_annotated.csv")
res <- fread(paste0(root, "interpretation/chromatin_analyses/results/results_df_with_promoterC_annotated.csv"))

#read eQTL enrichment from eQTLGen p2 to compare
res_eQTL_eqtlgen_ann_minimal <- fread(paste0(eqtldir, "output/Screen_results_minimal_annotated_eQTL_2025-12-23.tsv"))

#find eQTL files from eQTL Catalogue
datasets <- fread(paste0(eqtldir, "data/eQTL_catalogue/tabix_ftp_paths.tsv"))
eqtl_list <- list.files(paste0(eqtldir, "data/eQTL_catalogue/"), pattern = "QTD")

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
  # left_join(links_source, by = "link_id", relationship = 'many-to-many') %>% #this is where additional rows are created
  # mutate(source = case_when(source == "overlap" ~ "100kb window", #combine 100kb window and gene overlap into one category
  #                           source == "nonGWAS_eQTL" ~ "eQTL", #combine eQTLs into one category
  #                           str_detect(source, "GRN") ~ "GRN", #combine GRNs into one category
  #                           str_detect(source, "stream") ~ "up- and downstream", #combine up- and downstream into one category
  #                           is.na(source) ~ "not_prioritised",
  #                           .default = source),
  #        source = str_replace_all(source, "nonGWAS_", "")) %>%
  filter(pass_qc == TRUE) %>% #keep only results that pass QC and therefore were evaluated
  separate_rows(prioritization_sources, sep = ",\\s*") %>%
  distinct()

#make results into GRanges object 
gr.res <- res2 %>% 
  makeGRangesFromDataFrame(., seqnames.field="pert_chr", start.field="pert_start", end.field="pert_end", keep.extra.columns=T)
seqlevelsStyle(gr.res) <- "UCSC" #consistent gene formatting

#### COMPARE EQTL ENRICHMENT FOR EACH DATASET ####
all_log_odds <- list()
all_nr_eqtls <- list()

safe_log_odds_ratio <- function(a, b, c, d) {
  if (any(c(a, b, c, d) == 0)) {
    return(NA_real_)
  } else {
    return(log((a * d) / (b * c)))
  }
}

for (eqtl_file in eqtl_list) {
  
  #save dataset name
  dataset_id <- gsub(".credible_sets.tsv.gz", "", eqtl_file)
  dataset_info <- datasets %>% 
    filter(dataset_id == !!dataset_id) %>% 
    mutate(dataset_label = paste0(study_label, "_", sample_group)) %>%
    select(dataset_id, dataset_label)
  print(paste0("Now running eQTL enrichment for: ", dataset_info$dataset_label))
  
  #read eQTLs 
  eqtls <- fread(paste0(indir, eqtl_file))
  
  #keep only relevant columns
  eqtls <- eqtls %>%
    separate(variant, into = c("chr", "pos", "A1", "A2")) %>%
    select(gene_id, chr, pos, rsid, pip, pvalue, beta, se, z)
  
  #make eQTLs into GRanges object 
  gr.eqtls <- eqtls %>% 
    makeGRangesFromDataFrame(., seqnames.field="chr", start.field="pos", end.field="pos", keep.extra.columns=T)
  seqlevelsStyle(gr.eqtls) <- "UCSC"
  
  #overlap based on enhancer location and SNP location
  gr.res_overlap <- join_overlap_inner(gr.res, gr.eqtls, maxgap = 0)
  
  #keep only results where the eQTL gene also matches the result gene
  res_overlap_minimal <- gr.res_overlap %>%
    as.data.frame() %>%
    filter(gene_id.x == gene_id.y) %>%
    select(link_id, significant, rsid) %>%
    distinct() %>%
    group_by(link_id, significant) %>%
    summarise(variant = paste(unique(rsid), collapse = ","), .groups = "drop") #consolidate overlapping results to one line per peak-gene combination
  
  #annotate general results with eQTL overlaps, keep one line per enh-gene link
  res_eQTL_ann_minimal <- res2 %>%
    left_join(res_overlap_minimal, by = c("link_id", "significant")) %>%
    select(link_id, response_id, is_target, interaction_type, variant, significant) %>%
    group_by(link_id) %>%
    mutate(significant = any(significant)) %>%
    distinct()
  fwrite(res_eQTL_ann_minimal, paste0(outdir, "Screen_results_minimal_annotated_eQTL_", dataset_info$dataset_label, ".tsv"), sep = "\t", col.names = T, row.names = F)
  
  # NUMBER OF SIGNIFICANT EFFECTS OVERLAPPING AN EQTL
  nr_sign_eqtl_overlaps <- res_overlap_minimal %>%
    filter(significant == TRUE) %>%
    nrow()
  nr_eqtls <- data.frame(dataset_id = dataset_id, 
                         dataset_label = dataset_info$dataset_label,
                         nr_eQTLs = nr_sign_eqtl_overlaps)
  all_nr_eqtls[[eqtl_file]] <- nr_eqtls
  
  # ODDS RATIOS
  
  #significant vs. non-significant effects
  cont.table.sign <- table(factor(is.na(res_eQTL_ann_minimal$variant), levels = c(FALSE, TRUE)),
                           res_eQTL_ann_minimal$significant)
  
  sign_eqtl     <- cont.table.sign[1,2]
  nonsign_eqtl  <- cont.table.sign[1,1]
  sign_noeqtl   <- cont.table.sign[2,2]
  nonsign_noeqtl<- cont.table.sign[2,1]
  
  log_odds_ratio.sign <- safe_log_odds_ratio(sign_eqtl, nonsign_eqtl, sign_noeqtl, nonsign_noeqtl)
  
  #prioritised vs. non-prioritised effects
  cont.table.prio <- table(factor(is.na(res_eQTL_ann_minimal$variant), levels = c(FALSE, TRUE)),
                           res_eQTL_ann_minimal$is_target)
  
  sign_eqtl     <- cont.table.prio[1,2]
  nonsign_eqtl  <- cont.table.prio[1,1]
  sign_noeqtl   <- cont.table.prio[2,2]
  nonsign_noeqtl<- cont.table.prio[2,1]
  
  log_odds_ratio.prio <- safe_log_odds_ratio(sign_eqtl, nonsign_eqtl, sign_noeqtl, nonsign_noeqtl)
  
  #significant prioritised vs. non-significant prioritised effects
  cont.table.prio.sign <-  res_eQTL_ann_minimal %>%
    filter(is_target == TRUE) %>%
    mutate(any_variant = factor(!is.na(variant), levels = c(TRUE, FALSE))) %>%
    group_by(significant, any_variant) %>%
    tally() %>%
    ungroup() %>%
    complete(significant = unique(significant),
             any_variant = factor(c(TRUE, FALSE), levels = c(TRUE, FALSE)),
             fill = list(n = 0)) %>%
    mutate(any_variant = as.character(any_variant)) %>%
    pivot_wider(names_from = any_variant, values_from = n)
  
  sign_eqtl     <- cont.table.prio.sign[2, "TRUE", drop = TRUE]
  nonsign_eqtl  <- cont.table.prio.sign[1, "TRUE", drop = TRUE]
  sign_noeqtl   <- cont.table.prio.sign[2, "FALSE", drop = TRUE]
  nonsign_noeqtl<- cont.table.prio.sign[1, "FALSE", drop = TRUE]
  
  log_odds_ratio.prio.sign <- safe_log_odds_ratio(sign_eqtl, nonsign_eqtl, sign_noeqtl, nonsign_noeqtl)
  
  #Combine log odds ratios
  log_odds <- data.frame(dataset_id = dataset_id,
                         dataset_label = dataset_info$dataset_label,
                         log_OR_prioritised = log_odds_ratio.prio, 
                         log_OR_significant = log_odds_ratio.sign, 
                         log_OR_prioritised_significant = log_odds_ratio.prio.sign)

  all_log_odds[[eqtl_file]] <- log_odds
}

#Add the log odds and eQTL number datasets together
log_odds_df <- bind_rows(all_log_odds)
nr_eqtls_df <- bind_rows(all_nr_eqtls)
df <- log_odds_df %>%
  inner_join(nr_eqtls_df, by = c("dataset_id", "dataset_label")) %>%
  left_join(datasets, by = "dataset_id") %>%
  select(study_id, study_label, dataset_id, dataset_label, 
         log_OR_prioritised, log_OR_significant, log_OR_prioritised_significant, nr_eQTLs,
         sample_group, tissue_label, condition_label, sample_size) %>%
  mutate(simplified_tissue_label = gsub("\\s*\\([^\\)]*\\)", "", tissue_label),
         simplified_tissue_label = case_when(str_detect(simplified_tissue_label, "T.* cell") ~ "T cell",
                                             str_detect(simplified_tissue_label, "Treg") ~ "T cell",
                                             str_detect(simplified_tissue_label, "B.* cell") ~ "lymphoid (not T cell)",
                                             str_detect(simplified_tissue_label, "plasmablast") ~ "lymphoid (not T cell)",
                                             str_detect(simplified_tissue_label, "LCL") ~ "lymphoid (not T cell)",
                                             str_detect(simplified_tissue_label, "NK cell") ~ "lymphoid (not T cell)",
                                             str_detect(simplified_tissue_label, "plasma*") ~ "lymphoid (not T cell)",
                                             str_detect(simplified_tissue_label, "dendritic") ~ "myeloid",
                                             str_detect(simplified_tissue_label, "monocyte") ~ "myeloid",
                                             str_detect(simplified_tissue_label, "macrophage") ~ "myeloid",
                                             str_detect(simplified_tissue_label, "neur*") ~ "nervous system",
                                             str_detect(simplified_tissue_label, "microglia") ~ "nervous system",
                                             str_detect(simplified_tissue_label, "brain") ~ "nervous system",
                                             str_detect(simplified_tissue_label, "astrocyte") ~ "nervous system",
                                             str_detect(simplified_tissue_label, "neocortex") ~ "nervous system",
                                             str_detect(simplified_tissue_label, "ependymal") ~ "nervous system",
                                             str_detect(simplified_tissue_label, "floor plate") ~ "nervous system",
                                             str_detect(simplified_tissue_label, "nerve") ~ "nervous system",
                                             str_detect(simplified_tissue_label, "colon") ~ "digestive organs",
                                             str_detect(simplified_tissue_label, "intestine") ~ "digestive organs",
                                             str_detect(simplified_tissue_label, "stomach") ~ "digestive organs",
                                             str_detect(simplified_tissue_label, "esophagus") ~ "digestive organs",
                                             str_detect(simplified_tissue_label, "salivary") ~ "digestive organs",
                                             str_detect(simplified_tissue_label, "vagina") ~ "reproductive organs",
                                             str_detect(simplified_tissue_label, "uterus") ~ "reproductive organs",
                                             str_detect(simplified_tissue_label, "ovary") ~ "reproductive organs",
                                             str_detect(simplified_tissue_label, "placenta") ~ "reproductive organs",
                                             str_detect(simplified_tissue_label, "testis") ~ "reproductive organs",
                                             str_detect(simplified_tissue_label, "prostate") ~ "reproductive organs",
                                             str_detect(simplified_tissue_label, "breast") ~ "reproductive organs",
                                             str_detect(simplified_tissue_label, "liver") ~ "liver/kidney/spleen/pancreas",
                                             str_detect(simplified_tissue_label, "hepatocyte") ~ "liver/kidney/spleen/pancreas",
                                             str_detect(simplified_tissue_label, "kidney") ~ "liver/kidney/spleen/pancreas",
                                             str_detect(simplified_tissue_label, "spleen") ~ "liver/kidney/spleen/pancreas",
                                             str_detect(simplified_tissue_label, "pancrea*") ~ "liver/kidney/spleen/pancreas",
                                             str_detect(simplified_tissue_label, "heart") ~ "cardiovascular",
                                             str_detect(simplified_tissue_label, "artery") ~ "cardiovascular",
                                             str_detect(simplified_tissue_label, "synovium") ~ "skin/cartilage",
                                             str_detect(simplified_tissue_label, "cartilage") ~ "skin/cartilage",
                                             str_detect(simplified_tissue_label, "fibroblast") ~ "skin/cartilage",
                                             str_detect(simplified_tissue_label, "skin") ~ "skin/cartilage",
                                             str_detect(simplified_tissue_label, "adipose") ~ "adipose/muscle",
                                             str_detect(simplified_tissue_label, "muscle") ~ "adipose/muscle",
                                             str_detect(simplified_tissue_label, "thyroid") ~ "thyroid/adrenal/pituitary",
                                             str_detect(simplified_tissue_label, "adrenal") ~ "thyroid/adrenal/pituitary",
                                             str_detect(simplified_tissue_label, "pituitary") ~ "thyroid/adrenal/pituitary",
                                             str_detect(simplified_tissue_label, "hematopoietic") ~ "blood",
                                             str_detect(simplified_tissue_label, "platelet") ~ "blood",
                                             TRUE ~ simplified_tissue_label))
  
fwrite(df, paste0(outdir, "eQTL_enrichment_logOR_", Sys.Date(), ".tsv"), sep = "\t", col.names = T, row.names = F)
#df <- fread(paste0(outdir, "eQTL_enrichment_logOR_", Sys.Date(), ".tsv"))

##### add eQTLGen enrichment ##### 

#calculate the log(odds)
#significant vs. non-significant effects
cont.table.sign <- table(factor(res_eQTL_eqtlgen_ann_minimal$variant == "", levels = c(FALSE, TRUE)),
                         res_eQTL_eqtlgen_ann_minimal$significant)
log_odds_ratio.sign <- safe_log_odds_ratio(cont.table.sign[1,2], cont.table.sign[1,1], cont.table.sign[2,2], cont.table.sign[2,1])

#prioritised vs. non-prioritised effects
cont.table.prio <- table(factor(res_eQTL_eqtlgen_ann_minimal$variant == "", levels = c(FALSE, TRUE)),
                         res_eQTL_eqtlgen_ann_minimal$is_target)
log_odds_ratio.prio <- safe_log_odds_ratio(cont.table.prio[1,2], cont.table.prio[1,1], cont.table.prio[2,2], cont.table.prio[2,1])

#significant prioritised vs. non-significant prioritised effects
cont.table.prio.sign <-  res_eQTL_eqtlgen_ann_minimal %>%
  filter(is_target == TRUE) %>%
  mutate(any_variant = factor(variant != "", levels = c(TRUE, FALSE))) %>%
  group_by(significant, any_variant) %>%
  tally() %>%
  ungroup() %>%
  complete(significant = unique(significant), any_variant = factor(c(TRUE, FALSE), levels = c(TRUE, FALSE)), fill = list(n = 0)) %>%
  mutate(any_variant = as.character(any_variant)) %>%
  pivot_wider(names_from = any_variant, values_from = n)

log_odds_ratio.prio.sign <- safe_log_odds_ratio(cont.table.prio.sign[2, "TRUE", drop = TRUE], cont.table.prio.sign[1, "TRUE", drop = TRUE], cont.table.prio.sign[2, "FALSE", drop = TRUE], cont.table.prio.sign[1, "FALSE", drop = TRUE])

nr_eqtls <- res_eQTL_eqtlgen_ann_minimal %>%
  filter(significant == TRUE) %>%
  filter(variant != "") %>%
  nrow()
  
#Combine log odds ratios and nr of eQTLs
df_eqtlgen <- data.frame(study_id = "eQTLGen_p2",
                         study_label = "eQTLGen_p2",
                         dataset_id = "eQTLGen_p2_blood",
                         dataset_label = "eQTLGen_p2_blood",
                         log_OR_prioritised = log_odds_ratio.prio,
                         log_OR_significant = log_odds_ratio.sign, 
                         log_OR_prioritised_significant = log_odds_ratio.prio.sign,
                         nr_eQTLs = nr_eqtls,
                         sample_group = "blood",
                         tissue_label = "blood",
                         condition_label = "naive",
                         sample_size = 40000,
                         simplified_tissue_label = "blood")

df_ext <- rbind(df, df_eqtlgen)

#### PLOT ####

#colour palette to highlight T-cell and blood-related phenotypes
my_palette <- c(
  "T cell" = "firebrick1",
  "lymphoid (not T cell)" = "firebrick3",
  "myeloid" = "firebrick",
  "blood" = "firebrick4",
  "nervous system" = "goldenrod1",
  "adipose/muscle" = "darkgoldenrod1",
  "skin/cartilage" = "darkgoldenrod2",
  "thyroid/adrenal/pituitary" = "darkgoldenrod3",
  "cardiovascular" = "darkgoldenrod",
  "reproductive organs" = "darkgoldenrod4",
  "lung" = "burlywood2",
  "iPSC" = "burlywood",
  "digestive organs" = "burlywood3",
  "liver/kidney/spleen/pancreas" = "burlywood4"
)

#change order of factor for legend
df_ext <- df_ext %>%
  mutate(simplified_tissue_label = fct_relevel(df_ext$simplified_tissue_label, !!!names(my_palette)))
  
#Odds ratio of prioritised & significant results vs. total number of effects annotated with eQTLs
p <- ggplot(df_ext, aes(x = log_OR_prioritised_significant, y = nr_eQTLs)) +
  geom_vline(xintercept = 1, linetype = "dotted", colour = "grey") +
  geom_point(aes(colour = simplified_tissue_label)) +
  geom_text_repel(data = subset(df_ext, log_OR_prioritised_significant > 2.8 | nr_eQTLs > 45),
                  aes(label = tissue_label),
                  size = 2.5,
                  max.overlaps = Inf,
                  min.segment.length = 0,
                  segment.color = "grey50",
                  segment.size = 0.3) + 
  scale_colour_manual(values = my_palette) +
  theme_classic() +
  theme(
    legend.text = element_text(size = 7),      # Item text size
    legend.key.size = unit(0.4, "cm"),         # Box size
    legend.spacing.y = unit(0.4, "cm")         # Vertical spacing
  ) +
  labs(color = NULL) +
  ylab("Significant enhancer-gene pairs with eQTL\n") +
  xlab("\neQTL enrichment (log(OR)) among significant vs.\nnon-significant prioritised enhancer-gene pairs")
ggsave(plot = p, paste0(outdir, "eQTL_enrichment_across_datasets_", Sys.Date(), ".pdf"), 
       height = 5, width = 8)


