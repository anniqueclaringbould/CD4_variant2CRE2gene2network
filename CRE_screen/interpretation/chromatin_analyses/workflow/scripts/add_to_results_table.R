## add additional information to results table

# required packages
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(GenomicRanges)
  library(rtracklayer)
})

# Define functions ---------------------------------------------------------------------------------

# compute distance to target TSS or gene
compute_distance <- function(etp, feature = c("tss", "gene"),
                             fix = c("center", "none")) {
  
  # parse input arguments or set default value
  fix <- match.arg(fix)
  feature <- match.arg(feature)
  
  # GRanges object for all E-G pairs in ETP table using enhancers as coordinates
  pert <- makeGRangesFromDataFrame(etp, seqnames.field = "pert_chr", start.field = "pert_start",
                                   end = "pert_end", keep.extra.columns = TRUE)
  
  # resize perturbations to 1bp if specified
  if (fix != "none") {
    pert <- resize(pert, width = 1, fix = fix)
  }
  
  # GRanges object for all TSS or genes in ETP table
  if (feature == "tss") {
    feat <- makeGRangesFromDataFrame(etp, seqnames.field = "gene_chr", start.field = "gene_tss",
                                     end = "gene_tss", keep.extra.columns = TRUE)
  } else {
    feat <- makeGRangesFromDataFrame(etp, seqnames.field = "gene_chr", start.field = "gene_start",
                                     end = "gene_end", strand.field = "gene_strand",
                                     keep.extra.columns = TRUE)
  }
  
  # combine into a paired GRanges object containing enhancer and TSS/gene annotations for each pair
  eg_pairs <- Pairs(first = pert, second = feat)
  
  # compute distance between enhancers and provided TSS / gene annotations
  distance <- distance(eg_pairs)
  
  # get pairs where cCRE is upstream of the TSS depending on strand of the gene
  gene_strand <- as.vector(strand(second(eg_pairs)))
  upstream_pairs <- case_when(
    gene_strand == "+" & end(first(eg_pairs)) < start(second(eg_pairs)) ~ TRUE,
    gene_strand == "-" & start(first(eg_pairs)) > end(second(eg_pairs)) ~ TRUE,
    TRUE ~ FALSE
  )
  
  # set distance of upstream interactions to negative numbers
  distance[upstream_pairs] <- distance[upstream_pairs] * -1
  
  # add distance to etp table
  dist_col <- paste0("dist_to_", feature)
  etp <- etp %>% 
    mutate(!!dist_col := distance)
  
  return(etp)
  
}

# Add manual perturbation coordinates --------------------------------------------------------------

# load table with ETP results
etp <- read_csv(snakemake@input$etps, show_col_types = FALSE)

# remove redundant columns
etp <- select(etp, -c(enh_chr, enh_start, enh_end, seqid, strand, TSS_position))

# load manually annotated perturbation coordinates for promoter controls
pert_coords <- read_tsv(snakemake@input$pert_coords, show_col_types = FALSE)

# set manual coordinates for perturbations specified in manual coordinates file
etp_add_coords <- etp %>% 
  filter(grna_target %in% pert_coords$grna_target) %>% 
  select(-c(pert_chr, pert_start, pert_end)) %>% 
  left_join(pert_coords, by = "grna_target") %>% 
  select(all_of(c(colnames(etp), "ctrl_perturbation", "ctrl_type")))

# combine with E-G pairs that have coordinates
etp <- etp %>% 
  filter(!grna_target %in% pert_coords$grna_target) %>% 
  bind_rows(etp_add_coords) %>% 
  arrange(p_value)

# fill in control perturbation column
etp <- replace_na(etp, replace = list(ctrl_perturbation = FALSE))

# Add prioritization source ------------------------------------------------------------------------

# load processed E-G prioritization sources
prioritization <- read_csv(snakemake@input$eg_prioritization, show_col_types = FALSE)

# add unique identifiers for pairs to results
etp <- etp %>%
  mutate(grna_target_tmp = case_when(grna_target == "CD28" ~ "chr2:203706639-203706639", # add promoter controls as gene names (as in results)
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
  mutate(pair_uid = paste0(grna_target_tmp, "_", gene_id))

# add the prioritization sources to ETP table
etp <- etp %>%
  left_join(prioritization, by = "pair_uid") %>%
  replace_na(replace = list(prioritization_sources = "not prioritized")) %>% 
  select(-c(pair_uid, grna_target_tmp))

# Recompute distance to target gene ----------------------------------------------------------------

# load gene annotations
annot <- import(snakemake@input$gene_annot)

# calculate gene locus coordinates for each gene
genes <- split(annot, f = annot$gene_id)
genes <- unlist(range(genes))
genes$gene_id <- names(genes)

# replace gene coordinates with annotations used to recompute distance to genes
etp <- genes %>% 
  as.data.frame() %>% 
  select(gene_id, gene_chr = seqnames, gene_start = start, gene_end = end, gene_strand = strand) %>% 
  left_join(select(etp, -c(gene_chr, gene_start, gene_end, distance)), ., by = "gene_id")

# compute distance to target gene body (0 if in target gene)
etp <- compute_distance(etp, feature = "gene", fix = "center")

# Recompute distance to nearest target gene TSS ----------------------------------------------------

# calculate TSS coordinates for every annotated transcript
transcripts <- split(annot, f = annot$transcript_id)
tss <- promoters(unlist(range(transcripts)), upstream = 0, downstream = 1)
tss$transcript_id <- names(tss)

# add gene ids to tss list
transcript_gene_ids <- unique(mcols(annot)[, c("transcript_id", "gene_id")])
mcols(tss) <- merge(mcols(tss), transcript_gene_ids, by = "transcript_id", all.x = TRUE)

# replace TSS coordinates with annotations used to recompute distance to TSS
etp <- tss %>% 
  as.data.frame() %>% 
  select(gene_id, transcript_tss = transcript_id, gene_tss = start) %>% 
  left_join(select(etp, -c(gene_tss, TSS_dist)), ., by = "gene_id", relationship = "many-to-many")

# compute distance to all annotated TSSs per target gene
etp <- compute_distance(etp, feature = "tss", fix = "center")

# only retain closest TSS per perturbation - target gene pair
etp <- etp %>% 
  group_by(across(any_of(c("panel", "sample", "response_id", "grna_target")))) %>% 
  slice_min(abs(dist_to_tss), n = 1, with_ties = FALSE)

# reformat for output
etp <- etp %>% 
  select(-pval_scaled) %>% 
  relocate(response_id, grna_target, n_nonzero_trt, n_nonzero_cntrl, pass_qc, p_value,
           log_2_fold_change, significant, .before = 1) %>% 
  relocate(ctrl_perturbation, ctrl_type, target_genes, prioritization_sources,
           .before = "pert_chr") %>% 
  relocate(dist_to_gene, .before = "dist_to_tss") %>% 
  arrange(p_value, response_id, grna_target)
  
# save to output file
write_csv(etp, file = snakemake@output[[1]])
