##################################################################
## Script to analyse effects within TAD borders of CD4+ T cells ##
##################################################################

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

#directories
root <- "/g/steinmetz/project/otar_2063/enhancer_screen/"
indir <- paste0(root, "interpretation/hi-c_enrichment/data/")
outdir <- paste0(root, "interpretation/hi-c_enrichment/output/")

#### READ DATA ####

#read results, combined and annotated
#res <- fread("/g/steinmetz/gschwind/otar/manuscript_analyses/results/results_df_with_promoterC_annotated.csv")
res <- fread(paste0(root, "interpretation/chromatin_analyses/results/results_df_with_promoterC_annotated.csv"))

# for control perturbations which were in multiple panels, select the most significant one
res.ctrl <- res %>% 
  filter(ctrl_perturbation == TRUE) %>% 
  group_by(response_id, grna_target) %>% 
  slice_min(p_value, n = 1, with_ties = FALSE)

# add back to results table with discovery pairs
res <- res %>% 
  filter(ctrl_perturbation == FALSE) %>% 
  bind_rows(res.ctrl) %>% 
  arrange(desc(p_value)) 

#Filter the results to test
res.filt <- res %>%
  filter(pert_chr == gene_chr) %>% #keep enhancer-gene links on the same chromosome
  filter(pass_qc == TRUE) %>% #keep enhancer-gene links that were tested and pass QC
  filter(class_manuscript == "intergenic-distal" | true_enhancer == TRUE)

cell_types <- c("NS", "NU", "Th2U", "Th2S")

for (ct in cell_types){
  
  #read TAD boundaries
  tad <- fread(paste0(indir, "cons_", ct, ".bed"))
  #tad <- fread(paste0(indir, "cons_NS.bed")) #"cons_NU.bed" "cons_Th2U.bed" "cons_Th2S.bed" #to quickly check how the other TAD boundaries do
  
  #### ADJUST ####
  
  # Function to create regions within the TAD boundaries, for each chromosome
  create_regions <- function(df) {
    df %>%
      mutate(
        region_start = lag(end, default = 0),   # Use the previous boundary's end as the start of the region, default start at 0
        region_end = start                      # Current boundary's start is the region's end
      ) %>%
      select(chr, region_start, region_end) %>%
      mutate(region = paste0("region", row_number()))  # Create region names
  }
  
  # Apply the region creation separately for each chromosome
  tad_regions <- tad %>%
    select(chr = V1, start = V2, end = V3) %>%
    distinct() %>%
    group_by(chr) %>%
    do(create_regions(.)) %>%
    bind_rows( # Adjust the regions to handle the space after the last boundary in each chromosome
      tad %>%
        select(chr = V1, start = V2, end = V3) %>%
        group_by(chr) %>%
        summarise(
          region_start = max(end),  # Region starts at the last boundary's end
          region_end = Inf,         # Set as Inf for the last region
          region = paste0("region", n() + 1)  # Create the last region
        )
    ) %>%
    ungroup() %>%
    select(region_chr = chr, region_start, region_end, region)
  
  # Classify enhancers and targets into their respective regions
  classify_region <- function(chr, pos, tad_df) {
    tad_filtered <- tad_df %>% filter(region_chr == chr, region_start <= pos, region_end >= pos)
    if (nrow(tad_filtered) == 1) {
      return(tad_filtered$region)
    } else {
      return("boundary region")  # NA if no matching region it overlaps the boundary itself
    }
  }
  
  # Apply the region classification to the filtered dataset
  res.filt.ann <- res.filt %>%
    rowwise() %>%
    mutate(
      enhancer_region = classify_region(pert_chr, pert_start, tad_regions), #perturbation start is always upstream of perturbation end, so that's what we want to test
      target_region = classify_region(gene_chr, gene_start, tad_regions) #gene start is always upstream of gene end, so that's what we want to test
    ) %>%
    ungroup() %>%
    mutate(same_region = ifelse(enhancer_region == target_region, TRUE, FALSE)) # Check if enhancer and target are in the same region
  
  # Calculate the proportions of significant results within the same region vs across regions
  boundary_significance <- table(res.filt.ann$same_region, res.filt.ann$significant)
  
  #Fisher test
  fisher.test(boundary_significance)
  pval <- fisher.test(boundary_significance)$p.value
  odds_ratio <- fisher.test(boundary_significance)$estimate
  log_odds_ratio <- log(odds_ratio)
  
  #PLOT
  bin_levels <- c("TSS", "1500-10k", "10k-100k", "100k-1Mb", ">1Mb")
  
  # Bar plot of proportions while keeping non significant results
  res.filt.ann %>%
    filter(!is.na(same_region)) %>%
    filter(!is.na(distToTSS_bin)) %>%
    mutate(region = case_when(same_region == T ~ "within\nTADs",
                              same_region == F ~ "across\nTADs"),
           distToTSS_bin = fct_relevel(distToTSS_bin, bin_levels)) %>%
    ggplot(.) +
    geom_bar(aes(x = region, fill = significant), stat = "count", position = "fill", colour = "black") + 
    scale_fill_manual(values = c("grey70", "#0263a3")) +
    facet_grid(~distToTSS_bin) +
    labs(y = "Proportion of results\n", x = "") +
    theme_bw()
  ggsave(paste0(outdir, "Prop_results_TADs_by_distance_T_", ct, "_", Sys.Date(), ".pdf"), height = 4, width = 6)
  
  #Calculate proportions of significant results
  proportion_data <- res.filt.ann %>%
    group_by(same_region) %>%
    filter(!is.na(same_region)) %>%
    summarize(
      count_significant = sum(significant == TRUE),     # Count of significant results
      total_count = dplyr::n(),                                  # Total count of results
      pct_significant = count_significant / total_count * 100 # Calculate proportion
    ) %>%
    mutate(region = case_when(same_region == T ~ "within TADs",
                              same_region == F ~ "across TADs"))
  
  # Bar plot of proportions
  y_max <- max(proportion_data$pct_significant)
  
  ggplot(proportion_data, aes(x = region, y = pct_significant, fill = region, )) +
    geom_bar(stat = "identity", colour = "black") +
    scale_fill_manual(values = c("grey70", "#0263a3")) +
    geom_segment(
      aes(x = 1, xend = 2, y = y_max * 1.05, yend = y_max * 1.05),
      inherit.aes = FALSE) +
    annotate("text", x = 1.5, y = y_max * 1.12, label = paste0("OR = ", round(odds_ratio, digits = 1), "\np = ", signif(pval, 3))) +
    labs(y = "Percentage of results that is significant (%)\n") +
    theme_bw() +
    theme(legend.position = "none",
      axis.title.x = element_blank()) +
    scale_y_continuous(labels = scales::percent_format(scale = 1)) 
  
  ggsave(paste0(outdir, "Pct_significant_results_TADs_T_", ct, "_", Sys.Date(), ".pdf"), height = 5, width = 3.5)
  
  res.filt.ann %>%
    mutate(region = case_when(same_region == T ~ "within\nTADs",
                              same_region == F ~ "across\nTADs")) %>% 
    ggplot(.) +
    geom_bar(aes(x = region, fill = significant), stat = "count", position = "fill", colour = "black") + 
    scale_fill_manual(values = c("grey70", "#0263a3")) +
    geom_segment(
      aes(x = 1, xend = 2, y = 1.05, yend = 1.05),
      inherit.aes = FALSE) +
    annotate("text", x = 1.5, y = 1.12, label = paste0("OR = ", round(odds_ratio, digits = 1), "\np = ", signif(pval, 3))) +
    labs(y = "Proportion of results\n", x = "") +
    theme_bw()
  
  ggsave(paste0(outdir, "Prop_significant_results_TADs_T_", ct, "_", Sys.Date(), ".pdf"), height = 5, width = 3.5)
  
}

