## Annotate EPT results table with perturbation target type, regulatory interaction type and
## distance to TSS bin (for valid cis-acting pairs)

# save.image("RDA/annotate_etp_table.rda")
# stop()

# required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(GenomicRanges)
  library(rtracklayer)
})

# Define functions ---------------------------------------------------------------------------------

# overlap enhancer with promoters or gene bodies and get list of overlaps
overlap_enh_genes <- function(enh, features) {
  
  # overlap enhancers with features
  feature_ovl <- findOverlaps(enh, features)
  feature_enh <- tibble(
    enh = mcols(enh[queryHits(feature_ovl)])[["grna_target"]],
    id = mcols(features[subjectHits(feature_ovl)])[["id"]]
  )
  
  # make sure to retain unique enhancer-gene combinations
  feature_enh <- unique(feature_enh)
  
  return(feature_enh)
  
}

# Load and process input data ----------------------------------------------------------------------

# load table with ETP results
etp <- read_csv(snakemake@input$etps, show_col_types = FALSE)

# load gene annotations
genes <- import(snakemake@input$gene_annot)

# calculate gene locus coordinates for each gene
gene_loci <- split(genes, f = genes$gene_id)
gene_loci <- unlist(range(gene_loci))
gene_loci$id <- names(gene_loci)

# calculate all promoters for all annotated transcripts
transcripts <- split(genes, f = genes$transcript_id)
transcripts <- unlist(range(transcripts))
promoters <- promoters(transcripts, upstream = 1500, downstream = 1500)
promoters$id <- names(promoters)

# get gene ids for every transcript id
gene_id_per_transcript <- mcols(genes) %>% 
  as.data.frame() %>% 
  select(transcript_id, gene_id) %>% 
  distinct()

## Classify perturbations into promoter, genic and intergenic enhancers ----------------------------

# create GRanges object with enhancer perturbation coordinates
enh <- etp %>%
  filter(!is.na(pert_chr)) %>% 
  select(chr = pert_chr, start = pert_start, end = pert_end, grna_target) %>% 
  distinct() %>% 
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

# resize to center to be consistent with distance to TSS calculations
enh <- resize(enh, width = 1, fix = "center")

# overlap enhancers with gene promoters and gene bodies
enh_promoters <- overlap_enh_genes(enh, unique(promoters))
enh_genes <- overlap_enh_genes(enh, unique(gene_loci))

# get all enhancers overlapping promoters or gene bodies
prom_enh <- unique(enh_promoters$enh)
genic_enh <- setdiff(unique(enh_genes$enh), prom_enh)

# add information on promoter targeting and genic enhancers to etp results
etp <- etp %>% 
  mutate(element_type = case_when(
    grna_target %in% prom_enh ~ "promoter",
    grna_target %in% genic_enh ~ "genic",
    TRUE ~ "intergenic"
  ))

# add gene ids to promoters and genes overlapping enhancers
enh_promoters <- left_join(enh_promoters, gene_id_per_transcript, by = c("id" = "transcript_id"))

# add overlapping gene promoters to table
etp <- enh_promoters %>% 
  group_by(enh) %>% 
  summarize(overlapping_promoters = list(unique(gene_id))) %>% 
  left_join(etp, ., by = c("grna_target" = "enh"))

# add overlapping gene loci to table
etp <- enh_genes %>% 
  group_by(enh) %>% 
  summarize(overlapping_genes = list(unique(id))) %>% 
  left_join(etp, ., by = c("grna_target" = "enh"))

## Annotate ETP interaction types ------------------------------------------------------------------

# classify ETP interactions based on perturbed element and location of the target gene
etp <- etp %>% 
  rowwise() %>% 
  mutate(interaction_type = case_when(
    gene_id %in% overlapping_promoters | grna_target == response_id ~ "promoter-proximal",
    element_type == "promoter" & pert_chr == gene_chr ~ "promoter-distal",
    element_type == "promoter" ~ "promoter-trans",
    element_type == "genic" & dist_to_gene == 0 ~ "intragenic-proximal",
    element_type == "genic" & pert_chr == gene_chr ~ "intragenic-distal",
    element_type == "genic" ~ "intragenic-trans",
    element_type == "intergenic" & pert_chr == gene_chr ~ "intergenic-distal",
    element_type == "intergenic" ~ "intergenic-trans",
    TRUE ~ NA_character_
  )) %>% 
  ungroup()

## Bin ETP interactions into distance to TSS bins --------------------------------------------------

# annotate all valid cis-regulatory interactions based on their distance to TSS
etp <- etp %>% 
  mutate(distToTSS_bin = case_when(
    interaction_type == "promoter-proximal" ~ "TSS",
    (abs(dist_to_tss) > 1500 & abs(dist_to_tss) <= 1e4) ~ "1500-10k",
    (abs(dist_to_tss) > 1e4 & abs(dist_to_tss) <= 1e5) ~ "10k-100k",
    (abs(dist_to_tss) > 1e5 & abs(dist_to_tss) <= 1e6) ~ "100k-1Mb",
    abs(dist_to_tss) > 1e6 ~ ">1Mb",
    TRUE ~ NA_character_
  ))

## Label newly discovered distal cre-gene interactions ---------------------------------------------

# get novel enhancer-like CRE-gene interactions
etp_true_cis_interactions <- etp %>%
  filter(significant == TRUE)%>%
  filter(pert_chr == gene_chr) %>%
  filter(ctrl_perturbation == FALSE) %>%
  group_by(grna_target) %>%
  mutate(
    enhancer_like_interaction =
      case_when(
        # Rule 1: intragenic-distal + no intragenic-proximal in group
        interaction_type == "intragenic-distal" &
          !any(interaction_type == "intragenic-proximal") ~ TRUE,
        
        # Rule 2: promoter-distal + no promoter-proximal in group
        interaction_type == "promoter-distal" &
          !any(interaction_type == "promoter-proximal") ~ TRUE,
        
        # Rule 3: any intergenic-distal
        interaction_type == "intergenic-distal" ~ TRUE,
        
        # Otherwise FALSE
        TRUE ~ FALSE
      )
  )

# annotate each pair whether it's considered a new distal cis-regulatory interaction
etp <- etp_true_cis_interactions %>%
  ungroup() %>% 
  select(response_id, grna_target, enhancer_like_interaction) %>% 
  distinct() %>% 
  left_join(etp, ., by = c("response_id", "grna_target")) %>% 
  replace_na(replace = list(enhancer_like_interaction = FALSE))

## Reformat for output -----------------------------------------------------------------------------

# collapse list columns to strings and sort table according to significance
etp <- etp %>% 
  rowwise() %>% 
  mutate(overlapping_promoters = paste(overlapping_promoters, collapse = ","),
         overlapping_genes = paste(overlapping_genes, collapse = ",")) %>% 
  ungroup() %>% 
  arrange(p_value, response_id, grna_target)

# save annotated table to file
write_csv(etp, file = snakemake@output[[1]])
