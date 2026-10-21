## Apply one alternative filtering scenario to the ETP results and extract positive and negative
## CREs for chromatin feature comparisons. Re-implements the filters of annotate_etp_table.R and
## chromatin_features_etp_classes.Rmd with adjustable parameters (see config/alternative_filters.yml)

# save.image("RDA/alternative_filter.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(data.table)
  library(GenomicRanges)
  library(rtracklayer)
})

# get filter parameters for this scenario
filters <- snakemake@params$filters
promoter_window <- unlist(filters$promoter_window)
enhancer_like_classes <- unlist(filters$enhancer_like_classes)
proximal_rescue <- filters$proximal_rescue
filter_cres <- filters$filter_cres

# chromatin features to include in output
plot_features <- c("log_2_fold_change", "enh_signal_count_DNase-seq", "enh_signal_avg_H3K27ac",
                   "enh_signal_avg_H3K4me1", "enh_signal_avg_H3K4me3")

# Load and process input data ----------------------------------------------------------------------

# load table with ETP results (only required columns, as some columns contain very long strings)
etp_cols <- c("response_id", "grna_target", "pass_qc", "p_value", "significant", "ctrl_perturbation",
              "pert_chr", "pert_start", "pert_end", "gene_id", "gene_chr", "dist_to_gene",
              "dist_to_tss", plot_features)
etp <- fread(snakemake@input$etps, sep = ",", select = etp_cols)

# for control perturbations which were in multiple panels, select the most significant one
etp_ctrl <- etp %>%
  filter(ctrl_perturbation == TRUE) %>%
  group_by(response_id, grna_target) %>%
  slice_min(p_value, n = 1, with_ties = FALSE) %>%
  ungroup()

# add back to results table with discovery pairs and filter out any trans-chromosomal interactions
etp <- etp %>%
  filter(ctrl_perturbation == FALSE) %>%
  bind_rows(etp_ctrl) %>%
  filter(pert_chr == gene_chr)

# load gene annotations
genes <- import(snakemake@input$gene_annot)

# calculate gene locus coordinates for each gene
gene_loci <- unlist(range(split(genes, f = genes$gene_id)))

# calculate promoters for all annotated transcripts using the promoter window of this scenario
transcripts <- unlist(range(split(genes, f = genes$transcript_id)))
promoters <- promoters(transcripts, upstream = promoter_window[1],
                       downstream = promoter_window[2])

# get gene ids for every transcript id
gene_id_per_transcript <- mcols(genes) %>%
  as.data.frame() %>%
  select(transcript_id, gene_id) %>%
  distinct()

## Classify perturbations into promoter, genic and intergenic elements -----------------------------

# create GRanges object with element coordinates, resized to center to be consistent with distance
# to TSS calculations
enh <- etp %>%
  select(chr = pert_chr, start = pert_start, end = pert_end, grna_target) %>%
  distinct() %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
  resize(width = 1, fix = "center")

# overlap elements with promoters and gene bodies
prom_ovl <- findOverlaps(enh, promoters)
genes_ovl <- findOverlaps(enh, gene_loci)

# get all elements overlapping promoters or gene bodies
prom_enh <- unique(enh$grna_target[queryHits(prom_ovl)])
genic_enh <- setdiff(unique(enh$grna_target[queryHits(genes_ovl)]), prom_enh)

# get genes for each element that overlaps one of its promoters
enh_promoters <- tibble(grna_target = enh$grna_target[queryHits(prom_ovl)],
                        transcript_id = names(promoters)[subjectHits(prom_ovl)]) %>%
  left_join(gene_id_per_transcript, by = "transcript_id") %>%
  distinct(grna_target, gene_id)

# add element type to ETP results
etp <- etp %>%
  mutate(element_type = case_when(
    grna_target %in% prom_enh ~ "promoter",
    grna_target %in% genic_enh ~ "genic",
    TRUE ~ "intergenic"
  ))

## Annotate ETP interaction types and distance to TSS bins -----------------------------------------

# label pairs where the element overlaps a promoter of its target gene
etp <- etp %>%
  mutate(own_promoter = paste(grna_target, gene_id) %in%
           paste(enh_promoters$grna_target, enh_promoters$gene_id))

# classify cis-acting ETP interactions based on perturbed element and location of the target gene
etp <- etp %>%
  mutate(interaction_type = case_when(
    own_promoter | grna_target == response_id ~ "promoter-proximal",
    element_type == "promoter" ~ "promoter-distal",
    element_type == "genic" & dist_to_gene == 0 ~ "intragenic-proximal",
    element_type == "genic" ~ "intragenic-distal",
    TRUE ~ "intergenic-distal"
  ))

# distance to TSS bins, where the first bin starts at the upstream or downstream edge of the
# promoter window (e.g. "1500-10k" for symmetric or "2000/200-10k" for asymmetric windows)
first_bin <- paste0(paste(unique(promoter_window), collapse = "/"), "-10k")
bin_levels <- c("TSS", first_bin, "10k-100k", "100k-1Mb", ">1Mb")

# label pairs outside the promoter window based on strand-specific distance to TSS, where upstream
# distances are negative
etp <- etp %>%
  mutate(outside_window = dist_to_tss < -promoter_window[1] | dist_to_tss > promoter_window[2])

# annotate interactions based on their distance to TSS
etp <- etp %>%
  mutate(distToTSS_bin = case_when(
    interaction_type == "promoter-proximal" ~ bin_levels[1],
    (outside_window & abs(dist_to_tss) <= 1e4) ~ bin_levels[2],
    (abs(dist_to_tss) > 1e4 & abs(dist_to_tss) <= 1e5) ~ bin_levels[3],
    (abs(dist_to_tss) > 1e5 & abs(dist_to_tss) <= 1e6) ~ bin_levels[4],
    abs(dist_to_tss) > 1e6 ~ bin_levels[5],
    TRUE ~ NA_character_
  ))

## Label enhancer-like distal cre-gene interactions ------------------------------------------------

# for each element, check if it has significant or tested intragenic- or promoter-proximal pairs
proximal_pairs <- etp %>%
  filter(ctrl_perturbation == FALSE) %>%
  group_by(grna_target) %>%
  summarize(
    sig_intragenic = any(significant %in% TRUE & interaction_type == "intragenic-proximal"),
    tested_intragenic = any(pass_qc %in% TRUE & interaction_type == "intragenic-proximal"),
    sig_promoter = any(significant %in% TRUE & interaction_type == "promoter-proximal"),
    tested_promoter = any(pass_qc %in% TRUE & interaction_type == "promoter-proximal")
  )

# decide whether distal pairs of an element are allowed based on its proximal pairs
proximal_pairs <- proximal_pairs %>%
  mutate(
    allow_intragenic = case_when(
      proximal_rescue == "significant" ~ !sig_intragenic,
      proximal_rescue == "tested" ~ tested_intragenic & !sig_intragenic,
      proximal_rescue == "none" ~ TRUE
    ),
    allow_promoter = case_when(
      proximal_rescue == "significant" ~ !sig_promoter,
      proximal_rescue == "tested" ~ tested_promoter & !sig_promoter,
      proximal_rescue == "none" ~ TRUE
    )
  )

# label significant non-control pairs as enhancer-like based on interaction type (intragenic-
# proximal pairs are effects on host genes and therefore don't depend on proximal pairs)
etp <- etp %>%
  left_join(select(proximal_pairs, grna_target, allow_intragenic, allow_promoter),
            by = "grna_target") %>%
  replace_na(replace = list(allow_intragenic = FALSE, allow_promoter = FALSE)) %>%
  mutate(enhancer_like_interaction = significant %in% TRUE & ctrl_perturbation == FALSE & (
    interaction_type == "intergenic-distal" |
      (interaction_type == "intragenic-proximal" &
         "intragenic-proximal" %in% enhancer_like_classes) |
      (interaction_type == "intragenic-distal" & "intragenic-distal" %in% enhancer_like_classes &
         allow_intragenic) |
      (interaction_type == "promoter-distal" & "promoter-distal" %in% enhancer_like_classes &
         allow_promoter)
  ))

# annotate all rows of a pair as enhancer-like if any of them is
etp <- etp %>%
  group_by(response_id, grna_target) %>%
  mutate(enhancer_like_interaction = any(enhancer_like_interaction)) %>%
  ungroup()

## Extract positive and negative CREs --------------------------------------------------------------

# filter for pairs that passed QC
etp_filt <- filter(etp, !is.na(significant))

# filter for intergenic-distal, promoter-proximal and enhancer-like pairs
if (filter_cres == TRUE) {
  etp_filt <- etp_filt %>%
    filter(interaction_type %in% c("intergenic-distal", "promoter-proximal") |
             enhancer_like_interaction == TRUE)
}

# select any pairs that were detected as positives and categorize by distance
all_bins <- c(paste("hit", bin_levels), "NS")
positives <- etp_filt %>% 
  filter(significant == TRUE) %>% 
  mutate(category = paste("hit", distToTSS_bin)) %>%
  mutate(category = factor(category, levels = all_bins, ordered = TRUE) )

# select columns and features to plot
positives <- positives %>% 
  select(grna_target, significant, category, dist_to_tss, all_of(plot_features))

# if an enhancer has multiple positive hits, only pick the one closest to the TSS
positives <- positives %>% 
  group_by(grna_target) %>% 
  slice_min(abs(dist_to_tss), n = 1, with_ties = FALSE) %>% 
  select(-dist_to_tss)

# filter for elements that are not linked to any gene
negatives <- filter(etp_filt, significant == FALSE, !grna_target %in% positives$grna_target)

# only use intergenic-distal elements as negatives if CREs are filtered
if (filter_cres == TRUE) {
  negatives <- filter(negatives, interaction_type == "intergenic-distal")
}

# if an element is in multiple negative pairs, only pick the one closest to the TSS
negatives <- negatives %>%
  mutate(category = "NS") %>%
  select(grna_target, significant, category, dist_to_tss, all_of(plot_features)) %>% 
  group_by(grna_target) %>% 
  slice_min(abs(dist_to_tss), n = 1, with_ties = FALSE) %>% 
  select(-dist_to_tss)

# save all pairs with alternative filters to output file
write_tsv(etp, file = snakemake@output$pairs)

# save positives and negatives to output files
write_tsv(positives, file = snakemake@output$positives)
write_tsv(negatives, file = snakemake@output$negatives)
