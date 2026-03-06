###################################################
## Script to perform pathway enrichment analyses ##
###################################################

#Using https://rpubs.com/jrgonzalezISGlobal/enrichment as starting point

#### TODO 
#Databases: Cytopus Immune Sig (NatBiotech Peer), MSigDB, Reactome

#### SETUP ####
#libraries
library(tidyverse)
library(data.table)
library(ggplot2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(arrow)
library(rtracklayer)
library(BPCells)
library(msigdbr)
library(AnnotationDbi)

#directories
root <- "/g/steinmetz/project/otar_2063/"
indir <- paste0(root, "promoter_screen/interpretation/DE/output/combined/")
enrdir <- paste0(root, "promoter_screen/interpretation/enrichment/")
mofadir <- paste0(root, "promoter_screen/interpretation/MOFA/")
outdir <- paste0(enrdir, "output/")

select <- dplyr::select
first <- dplyr::first


#### READ DATA ####
#promoter screen results
#res <- fread("/g/steinmetz/project/otar_2063/promoter_screen/interpretation/diff_expr_analyses/process_trans_analysis/gene_level_empirical_pvalues/significant_trans_pairs_per_gene_fdr.tsv.gz")
#res <- fread("/g/steinmetz/project/otar_2063/promoter_screen/interpretation/diff_expr_analyses/process_trans_analysis/gene_level_empirical_pvalues/upsample_nt_pvals/discovery_pairs_empirical_pvalues.tsv.gz")
#res <- fread(paste0(indir, "DEG_single_cell_limma_ebayes_split_NT_removed_all_covariates_prefiltered_genes_gene_combined_p01_2026-01-21.tsv.gz"))
#res <- fread(paste0(indir, "DEG_single_cell_limma_ebayes_split_NT_removed_all_covariates_0.1pct_per_split_gene_combined_p01_2026-01-20.tsv.gz"))

#pbulk_deseq_files <- list.files("/g/stegle/schrod/code/TCell/pseudobulk_deseq2_chunks50/", full.names = TRUE)
pbulk_deseq_files <- list.files("/g/stegle/schrod/code/TCell/pseudobulk_deseq2_chunks50_min3bulks_ribosomal/", full.names = TRUE)
dm <- rbindlist(lapply(pbulk_deseq_files, fread))
pval_threshold <- 0.1
res <- dm %>%
  transmute(
    perturbation = contrast,
    target       = variable,
    pvalue       = p_value,
    adj_pvalue   = adj_p_value,
    significant  = adj_p_value < pval_threshold,
    is_target    = perturbation == target,
    logFC       = log_fc
  )

#enhancer results annotated by the traits
enh <- fread(paste0(enrdir, "input/cre_perturbation_screen_results_traits.csv"))

#MOFA 
#mofa_pert <- fread(paste0(mofadir, "T_cell_project/MOFA_Leiden_clusters.csv"))
mofa_pert <- fread(paste0(mofadir, "MOFA_Jana/MOFA_Leiden_clusters_perturbations.csv"))
mofa_ds <- fread(paste0(mofadir, "MOFA_Jana/MOFA_Leiden_clusters_downstream_genes.csv"))

#open targets database
ot <- open_dataset(paste0(enrdir, "input/association_by_datasource_direct")) |> collect()
ot_diseases <- open_dataset(paste0(enrdir, "input/disease")) |> collect()

#GTF
gtf <- readGFF(paste0(root, "/references/alignment_reference/hg38/gencode_gtf/default/hg38.gtf"))


#### PREPARATION ####

#CONVERSION FUNCTIONS
# symbols_to_entrez <- function(genes) {
#   unique(unlist(
#     mget(genes, envir = org.Hs.egALIAS2EG, ifnotfound = NA)
#   ))
# }
# 
# symbols_to_entrez <- function(genes) {
#   mapIds(
#     org.Hs.eg.db,
#     keys = genes,
#     keytype = "SYMBOL",
#     column = "ENTREZID",
#     multiVals = "first"
#   ) |> na.omit() |> unique()
# }

symbols_to_entrez <- function(genes) {
  
  genes <- unique(genes)
  
  # ---- 1. Strict SYMBOL mapping only for valid symbols ----
  valid_symbols <- intersect(genes, keys(org.Hs.eg.db, keytype = "SYMBOL"))
  
  sym_map <- rep(NA_character_, length(genes))
  names(sym_map) <- genes
  
  if (length(valid_symbols) > 0) {
    sym_map[valid_symbols] <- mapIds(
      org.Hs.eg.db,
      keys = valid_symbols,
      keytype = "SYMBOL",
      column = "ENTREZID",
      multiVals = "first"
    )
  }
  
  # ---- 2. Fallback to ALIAS for remaining genes ----
  remaining <- names(sym_map)[is.na(sym_map)]
  
  if (length(remaining) > 0) {
    
    valid_alias <- intersect(remaining, keys(org.Hs.eg.db, keytype = "ALIAS"))
    
    if (length(valid_alias) > 0) {
      sym_map[valid_alias] <- mapIds(
        org.Hs.eg.db,
        keys = valid_alias,
        keytype = "ALIAS",
        column = "ENTREZID",
        multiVals = "first"
      )
    }
  }
  
  unique(na.omit(sym_map))
}

symbols_to_ensembl <- function(genes, gtf) {
  as.data.frame(genes) %>%
    setNames("gene_name") %>%
    left_join(gtf, by = "gene_name") %>%
    select(gene_id) %>%
    distinct() %>%
    pull()
}

ensembl_to_symbols <- function(genes, gtf) {
  as.data.frame(genes) %>%
    setNames("gene_id") %>%
    left_join(gtf, by = "gene_id") %>%
    select(gene_name) %>%
    #distinct() %>% #to use it in a dataframe where the same genes occur more than once
    pull()
}

#GENES

#Adjust GTF to keep only ENSG IDs and gene names
gtf <- gtf %>%
  mutate(gene_id = gsub("\\..*", "", gene_id)) %>%
  mutate(gene_id = gsub("_.*", "", gene_id)) %>%
  dplyr::select(gene_id, gene_name) %>%
  distinct()

#List of all tested genes ("gene universe")
#all_genes_symbol <- unique(res$Target) #Limma results
all_genes_symbol <- unique(res$target) #deseq results
# all_genes_symbol <- res %>% #Gene-level p-value adjusted SCEPTRE results
#   select(response_id) %>%
#   distinct() %>%
#   mutate(response_id = ensembl_to_symbols(gsub("\\..*", "", response_id), gtf)) %>%
#   pull(response_id)

#Convert to ENTREZ and ENSG IDs for enrichments
all_genes_entrez <- symbols_to_entrez(all_genes_symbol)
all_genes_ensembl <- symbols_to_ensembl(all_genes_symbol, gtf) #Limma and deseq results
# all_genes_ensembl <- res %>% #Gene-level p-value adjusted SCEPTRE results
#   select(response_id) %>%
#   distinct() %>%
#   mutate(response_id = gsub("\\..*", "", response_id)) %>%
#   pull(response_id)

#PREP OPEN TARGETS DATABASE

#keep only disease ID and description
ot_diseases <- ot_diseases %>%
  select(id, name, description)

#Prep disease to gene and name conversions
disease2name <- ot_diseases[, c("id", "name")]
disease2gene <- ot[, c("diseaseId", "targetId")]

#PREP IMMUNESIGDB DATABASE

#Load Immune Signatures from MSigDB
immunesigdb <- msigdbr(species = "Homo sapiens",
                        category = "C7",
                        subcategory = "IMMUNESIGDB")

#Prep disease to gene and name conversions
gs2name <- immunesigdb[, c("gs_id", "gs_description")]
gs2gene <- immunesigdb[, c("gs_id", "ensembl_gene")]

#PREP HALLMARK DATABASE

#Load Hallmark gene sets from MSigDB
hmdb <- msigdbr(species = "Homo sapiens",
                category = "H")

#Prep disease to gene and name conversions
gs2name.hm <- hmdb[, c("gs_id", "gs_description")]
gs2gene.hm <- hmdb[, c("gs_id", "ensembl_gene")]

#ENRICHMENT FUNCTIONS

#GO enrichment
run_go_enrichment <- function(fg_genes_entrez, bg_genes_entrez,
                              ont = "ALL", pcut = 1, padj = "BH") {
  enr <- enrichGO(
    gene = fg_genes_entrez,
    universe = bg_genes_entrez,
    ont = ont,
    OrgDb = "org.Hs.eg.db",
    keyType = "ENTREZID",
    readable = TRUE,
    pvalueCutoff = pcut,
    pAdjustMethod = padj,
    pool = FALSE
  )
  
  # enr_simple <- simplify(enr, cutoff = 0.7,
  #                        by = "p.adjust",
  #                        select_fun = min)
}

#KEGG enrichment
run_kegg_enrichment <- function(fg_genes_entrez, bg_genes_entrez,
                                organism = "hsa", pcut = 1, padj = "BH") {
  enr <- enrichKEGG(
    gene = fg_genes_entrez,
    universe = bg_genes_entrez,
    organism = organism,
    keyType = "kegg",
    pvalueCutoff = pcut,
    pAdjustMethod = padj
  )
}

#Open Targets enrichment
run_ot_enrichment <- function(fg_genes_ensembl, bg_genes_ensembl,
                              disease2gene, disease2name,
                              pcut = 1, padj = "BH") {
  enr <- enricher(
    gene = fg_genes_ensembl,
    universe = bg_genes_ensembl,
    TERM2GENE = disease2gene,
    TERM2NAME = disease2name,
    pvalueCutoff = pcut,
    pAdjustMethod = padj
  )
}

#ImmuneSigDB enrichment
run_immunesigdb_enrichment <- function(fg_genes_ensembl, bg_genes_ensembl,
                                       gs2gene, gs2name,
                                       pcut = 1, padj = "BH") {
  enr <- enricher(
    gene = fg_genes_ensembl,
    universe = bg_genes_ensembl,
    TERM2GENE = gs2gene,
    TERM2NAME = gs2name,
    pvalueCutoff = pcut,
    pAdjustMethod = padj
  )
}

#Hallmark geneset enrichment
run_hmdb_enrichment <- function(fg_genes_ensembl, bg_genes_ensembl,
                                gs2gene.hm, gs2name.hm,
                                pcut = 1, padj = "BH") {
  enr <- enricher(
    gene = fg_genes_ensembl,
    universe = bg_genes_ensembl,
    TERM2GENE = gs2gene.hm,
    TERM2NAME = gs2name.hm,
    pvalueCutoff = pcut,
    pAdjustMethod = padj
  )
}

#FUNCTIONS TO RUN ALL ENRICHMENTS
run_all_enrichments <- function(
    genes_symbol,
    gtf,
    all_genes_entrez,
    all_genes_ensembl,
    disease2gene,
    disease2name,
    outdir,
    label
) {
  
  if (length(genes_symbol) == 0) return(NULL)
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  #Save number of genes
  n_genes <- length(unique(genes_symbol))
  plot_title <- paste0(label, " (n = ", n_genes, " genes)")
  
  #Convert genes to entrez or ensembl IDs
  fg_entrez  <- symbols_to_entrez(genes_symbol)
  fg_ensembl <- symbols_to_ensembl(genes_symbol, gtf)
  
  ## ---- GO ----
  message("GO enrichment")

  try({
    go <- run_go_enrichment(fg_entrez, all_genes_entrez)
    go.df <- as.data.frame(go)

    if (nrow(go.df) > 0) {
      fwrite(go.df, file.path(outdir, paste0("GO_", Sys.Date(), ".tsv")), sep = "\t")

      ggsave(
        file.path(outdir, paste0("GO_", Sys.Date(), ".pdf")),
        dotplot(go, showCategory = 10, color = "p.adjust",
                title = plot_title, label_format = 50),
        width = 8, height = 6
      )
    } else {
      message("No enrichments for GO")
    }
  }, silent = TRUE)

  # ## ---- KEGG ----
  # message("KEGG enrichment")
  # 
  # try({
  #   kegg <- run_kegg_enrichment(fg_entrez, all_genes_entrez)
  #   kegg.df <- as.data.frame(kegg)
  # 
  #   if (nrow(kegg.df) > 0) {
  #     fwrite(kegg.df, file.path(outdir, paste0("KEGG_", Sys.Date(), ".tsv")), sep = "\t")
  # 
  #     ggsave(
  #       file.path(outdir, paste0("KEGG_", Sys.Date(), ".pdf")),
  #       dotplot(kegg, showCategory = 10, color = "p.adjust",
  #               title = plot_title, label_format = 50),
  #       width = 8, height = 6
  #     )
  #   } else {
  #     message("No enrichments for KEGG")
  #   }
  # }, silent = TRUE)
  # 
  # ## ---- Open Targets ----
  # message("Open Targets enrichment")
  # 
  # try({
  #   ot <- run_ot_enrichment(
  #     fg_ensembl, all_genes_ensembl,
  #     disease2gene, disease2name
  #   )
  #   ot.df <- as.data.frame(ot)
  # 
  #   if (nrow(ot.df) > 0) {
  #     fwrite(ot.df, file.path(outdir, paste0("OpenTargets_", Sys.Date(), ".tsv")), sep = "\t")
  # 
  #     ggsave(
  #       file.path(outdir, paste0("OpenTargets_", Sys.Date(), ".pdf")),
  #       dotplot(ot, showCategory = 10, color = "p.adjust",
  #               title = plot_title, label_format = 50),
  #       width = 8, height = 6
  #     )
  #   } else {
  #     message("No enrichments for Open Targets")
  #   }
  # }, silent = TRUE)
  # 
  # ## ---- ImmuneSigDB ----
  # message("ImmuneSigDB enrichment")
  # 
  # try({
  #   imsig <- run_immunesigdb_enrichment(
  #     fg_ensembl, all_genes_ensembl,
  #     gs2gene, gs2name
  #   )
  #   imsig.df <- as.data.frame(imsig)
  # 
  #   if (nrow(imsig.df) > 0) {
  #     fwrite(imsig.df, file.path(outdir, paste0("ImmuneSigDB_", Sys.Date(), ".tsv")), sep = "\t")
  # 
  #     ggsave(
  #       file.path(outdir, paste0("ImmuneSigDB_", Sys.Date(), ".pdf")),
  #       dotplot(imsig, showCategory = 10, color = "p.adjust",
  #               title = plot_title, label_format = 50) +
  #         theme(axis.text.y = element_text(size = 8)),
  #       width = 8, height = 6
  #     )
  #   } else {
  #     message("No enrichments for ImmuneSigDB")
  #   }
  # }, silent = TRUE)
  # 
  # ## ---- Hallmark ----
  # message("Hallmark gene set enrichment")
  # 
  # try({
  #   hm <- run_hmdb_enrichment(
  #     fg_ensembl, all_genes_ensembl,
  #     gs2gene.hm, gs2name.hm
  #   )
  #   hm.df <- as.data.frame(hm)
  #   
  #   if (nrow(hm.df) > 0) {
  #     fwrite(hm.df, file.path(outdir, paste0("Hallmark_", Sys.Date(), ".tsv")), sep = "\t")
  #     
  #     ggsave(
  #       file.path(outdir, paste0("Hallmark_", Sys.Date(), ".pdf")),
  #       dotplot(hm, showCategory = 10, color = "p.adjust",
  #               title = plot_title, label_format = 50) +
  #         theme(axis.text.y = element_text(size = 8)),
  #       width = 8, height = 6
  #     )
  #   } else {
  #     message("No enrichments for Hallmark")
  #   }
  # }, silent = TRUE)
  
}

run_enrichment_sets <- function(
    gene_sets,
    prefix,
    ...
) {
  
  for (name in names(gene_sets)) {
    
    message("Running enrichment for ", name)
    
    out <- file.path(prefix, paste0(name, "_downstream_genes"))
    
    run_all_enrichments(
      genes_symbol = gene_sets[[name]],
      outdir = out,
      label = name,
      ...
    )
  }
}


#### ENRICHMENTS ####

#### SET OF GENES LINKED TO INDIVIDUAL PERTURBATIONS ####

# #List of genes downstream of perturbations that affect multiple genes

# #Limma results
# genes_by_perturbation <- res %>%
#   filter(significant) %>%
#   group_by(Perturbation) %>%
#   summarise(genes = list(Target), .groups = "drop") %>%
#   filter(lengths(genes) >= 800) %>%
#   deframe()

# #Gene-level p-value adjusted SCEPTRE results
# genes_by_perturbation <- res %>%
#   filter(pval_adj_per_gene < 0.1) %>%
#   mutate(response_id = ensembl_to_symbols(gsub("\\..*", "", response_id), gtf)) %>%
#   group_by(grna_target) %>%
#   summarise(genes = list(response_id), .groups = "drop") %>%
#   filter(lengths(genes) >= 100) %>%
#   deframe()

#Deseq results
genes_by_perturbation <- res %>%
  filter(significant) %>%
  group_by(perturbation) %>%
  summarise(genes = list(target), .groups = "drop") %>%
  filter(lengths(genes) >= 500) %>%
  deframe()

#Run enrichments for all perturbations with many targets
run_enrichment_sets(
  genes_by_perturbation,
  prefix = paste0(outdir, "perturbations"),
  gtf,
  all_genes_entrez,
  all_genes_ensembl,
  disease2gene,
  disease2name
)

#EXAMPLE LOCI

#select promoter screen results that are significant and with abs(logFC) > 0.2 for hop2
res.filt <- res %>%
  filter(significant == TRUE) %>%
  filter(abs(logFC) > 0.2)

#combine with results again to get hop2 genes
res_hop2 <- res %>%
  filter(significant == TRUE) %>%
  select(perturbation, target_hop1 = target, is_target_hop1 = is_target) %>%
  left_join(res.filt, by = c("target_hop1" = "perturbation"), relationship = "many-to-many") %>%
  select(perturbation, target_hop1, is_target_hop1, target_hop2 = target, is_target_hop2 = is_target) %>%
  distinct()

#Function to extract genes for one perturbation and run enrichment
run_enrichment_perturbation <- function(gene_name) {
  
  genes_vec <- res %>%
    filter(significant == TRUE) %>%
    group_by(perturbation) %>%
    filter(perturbation == gene_name) %>%
    summarise(
      genes = list(target),
      .groups = "drop"
    ) %>%
    deframe()
  
  if (length(genes_vec) == 0) {
    message("No genes found for ", gene_name)
    return(NULL)
  }
  
  run_enrichment_sets(
    genes_vec,
    prefix = paste0(outdir, "perturbations"),
    gtf,
    all_genes_entrez,
    all_genes_ensembl,
    disease2gene,
    disease2name
  )
}

#Function to extract genes for one perturbation and run enrichment on all hop 1 + hop 2 genes
run_enrichment_perturbation_hop2 <- function(gene_name) {
  
  genes_vec <- res_hop2 %>%
    group_by(perturbation) %>%
    filter(perturbation == gene_name) %>%
    summarise(
      genes = list(sort(unique(c(target_hop1, target_hop2)))),
      .groups = "drop"
    ) %>%
    deframe()
  
  if (length(genes_vec) == 0) {
    message("No genes found for ", gene_name)
    return(NULL)
  }
  
  run_enrichment_sets(
    genes_vec,
    prefix = paste0(outdir, "perturbations_hop2"),
    gtf,
    all_genes_entrez,
    all_genes_ensembl,
    disease2gene,
    disease2name
  )
}

#DEXI
run_enrichment_perturbation("DEXI")
run_enrichment_perturbation_hop2("DEXI")

#DEXI hub genes from hop1
hop1_hub_genes_DEXI <- res_hop2 %>%
  group_by(perturbation) %>%
  filter(perturbation == "DEXI") %>%
  group_by(target_hop1) %>% 
  tally() %>%
  arrange(-n) %>%
  filter(n > 5)

for (g in hop1_hub_genes_DEXI$target_hop1) {
  run_enrichment_perturbation(g)
}

#TYK2
run_enrichment_perturbation("TYK2")
run_enrichment_perturbation_hop2("TYK2")

#TYK2 hub genes from hop1
hop1_hub_genes_TYK2 <- res_hop2 %>%
  group_by(perturbation) %>%
  filter(perturbation == "TYK2") %>%
  group_by(target_hop1) %>% 
  tally() %>%
  arrange(-n) %>%
  filter(n > 5)

for (g in hop1_hub_genes_TYK2$target_hop1) {
  run_enrichment_perturbation(g)
}

#### SET OF GENES LINKED TO ENHANCERS FROM ONE DISEASE #### 

# #Limma results
# # #Combine enhancer and promoter results
# # comb <- enh %>%
# #   filter(significant == T) %>% #keep only significant enhancer effects
# #   select(grna_target, response_id, is_target_enh = is_target, p_value, log_2_fold_change, ctrl_perturbation, ctrl_type, element_type, interaction_type, enhancer_like_interaction, trait) %>%
# #   inner_join(res, by = c("response_id" = "Perturbation"), relationship = "many-to-many") %>%
# #   filter(significant == T) %>% #keep only significant promoter effects
# #   select(grna_target, target_enh = response_id, is_target_enh, target_prom = Target, is_target_prom = is_target, 
# #          p_value_enh = p_value, log_2_fold_change_enh = log_2_fold_change, 
# #          ctrl_perturbation, ctrl_type, element_type, interaction_type, enhancer_like_interaction,
# #          p_value_prom = P.Value, p_value_adj_prom = adj.P.Val, log_2_fold_change_prom = logFC, trait) %>%
# #   distinct()
# # 
# # #Adjust to select relevant columns and group by trait (disease)
# # comb_info <- comb %>%
# #   select(grna_target, target_enh, is_target_enh, target_prom, is_target_prom, ctrl_perturbation, ctrl_type, trait) %>%
# #   distinct() %>%
# #   separate_rows(trait, sep = ",") %>%
# #   group_by(trait) %>%
# #   add_tally()
# # 
# # table(is.na(comb_info$trait), comb_info$ctrl_perturbation)
# # #shows that all values where trait is NA, it's because the 'enhancer' was actually a control element
# # 
# # #List of diseases to use for filtering of downstream results
# # genes_by_disease <- comb_info %>%
# #   filter(!is.na(trait)) %>%
# #   group_by(trait) %>%
# #   summarise(genes = list(target_prom), .groups = "drop") %>%
# #   deframe()
# 
# #Gene-level p-value adjusted SCEPTRE results
# #Combine enhancer and promoter results
# comb <- enh %>%
#   filter(significant == T) %>% #keep only significant enhancer effects
#   select(grna_target, response_id, is_target_enh = is_target, p_value, log_2_fold_change, ctrl_perturbation, ctrl_type, element_type, interaction_type, enhancer_like_interaction, trait) %>%
#   inner_join(res, by = c("response_id" = "grna_target"), relationship = "many-to-many") %>%
#   filter(pval_adj_per_gene < 0.1) %>% #keep only significant promoter effects
#   mutate(target_prom = ensembl_to_symbols(gsub("\\..*", "", response_id.y), gtf)) %>%
#   select(grna_target, target_enh = response_id, is_target_enh, target_prom, 
#          p_value_enh = p_value.x, log_2_fold_change_enh = log_2_fold_change.x, 
#          ctrl_perturbation, ctrl_type, element_type, interaction_type, enhancer_like_interaction,
#          p_value_prom = p_value.y, p_value_adj_prom = pval_adj_per_gene, log_2_fold_change_prom = log_2_fold_change.y, trait) %>%
#   distinct()
# 
# #Adjust to select relevant columns and group by trait (disease)
# comb_info <- comb %>%
#   select(grna_target, target_enh, is_target_enh, target_prom, ctrl_perturbation, ctrl_type, trait) %>%
#   distinct() %>%
#   separate_rows(trait, sep = ",") %>%
#   group_by(trait) %>%
#   add_tally()
# 
# table(is.na(comb_info$trait), comb_info$ctrl_perturbation)
# #shows that all values where trait is NA, it's because the 'enhancer' was actually a control element
# 
# #List of diseases to use for filtering of downstream results
# genes_by_disease <- comb_info %>%
#   filter(!is.na(trait)) %>%
#   group_by(trait) %>%
#   summarise(genes = list(target_prom), .groups = "drop") %>%
#   deframe()


#Deseq results
#Combine enhancer and promoter results
comb <- enh %>%
  filter(significant == T) %>% #keep only significant enhancer effects
  #filter(gene_chr == pert_chr) %>% #NEVER RAN THIS YET BUT WOULD BE BETTER
  select(grna_target, response_id, is_target_enh = is_target, p_value, log_2_fold_change, 
         ctrl_perturbation, ctrl_type, element_type, interaction_type, enhancer_like_interaction, 
         distToTSS_bin, trait) %>%
  inner_join(res, by = c("response_id" = "perturbation"), relationship = "many-to-many") %>%
  filter(significant == T) %>% #keep only significant promoter effects
  select(grna_target, target_enh = response_id, is_target_enh, target_prom = target, is_target_prom = is_target, 
         p_value_enh = p_value, log_2_fold_change_enh = log_2_fold_change,
         ctrl_perturbation, ctrl_type, element_type, interaction_type, enhancer_like_interaction, distToTSS_bin,
         p_value_prom = pvalue, p_value_adj_prom = adj_pvalue, log_2_fold_change_prom = logFC, trait) %>%
  distinct()

#Adjust to select relevant columns and group by trait (disease)
comb_info <- comb %>%
  select(grna_target, target_enh, is_target_enh, target_prom, is_target_prom, ctrl_perturbation, ctrl_type, trait) %>%
  distinct() %>%
  separate_rows(trait, sep = ",") %>%
  mutate(trait = case_when(trait == "Primary_biliary_cirrhosis" ~ "Primary_biliary_cholangitis",
                           TRUE ~ trait)) %>% #cirrhosis older name for same disease, rename and combine
  group_by(trait) %>%
  add_tally()

table(is.na(comb_info$trait), comb_info$ctrl_perturbation)
#shows that all values where trait is NA, it's because the 'enhancer' was actually a control element

#List of diseases to use for filtering of downstream results
genes_by_disease <- comb_info %>%
  filter(!is.na(trait)) %>%
  group_by(trait) %>%
  summarise(genes = list(target_prom), .groups = "drop") %>%
  mutate(n_genes = lengths(genes)) %>%
  filter(n_genes > 3) %>% #remove diseases where the gene list only contains 3 genes or less
  select(-n_genes) %>%
  deframe()

#Run enrichments
run_enrichment_sets(
  genes_by_disease,
  prefix = paste0(outdir, "diseases"),
  gtf,
  all_genes_entrez,
  all_genes_ensembl,
  disease2gene,
  disease2name
)

#### SET OF GENES + DOWNSTREAM GENES ("HOP 2 GENES") LINKED TO ENHANCERS FROM ONE DISEASE #### 

#select promoter screen results that are significant and with abs(logFC) > 0.2
res.filt <- res %>%
  filter(significant == TRUE) %>%
  filter(abs(logFC) > 0.2)

#Adjust list of enhancers and link to downstream genes (hop 1) + downstream genes (hop 2)
genes_by_disease_hop2 <- comb_info %>%
  filter(!is.na(trait)) %>%
  left_join(res.filt,
            by = c("target_prom" = "perturbation"),
            relationship = "many-to-many") %>%
  select(trait, target_prom, target_prom_hop2 = target) %>%
  pivot_longer(cols = c(target_prom, target_prom_hop2),
               values_to = "gene") %>%
  filter(!is.na(gene)) %>%
  distinct(trait, gene) %>%
  group_by(trait) %>%
  summarise(genes = list(sort(gene)), .groups = "drop") %>%
  mutate(n_genes = lengths(genes)) %>%
  filter(n_genes > 3) %>%
  select(-n_genes) %>%
  deframe()

#Run enrichments
run_enrichment_sets(
  genes_by_disease_hop2,
  prefix = paste0(outdir, "diseases_hop2"),
  gtf,
  all_genes_entrez,
  all_genes_ensembl,
  disease2gene,
  disease2name
)


#### SET OF GENES LINKED TO ONE ENHANCER #### 

#Adjust to select relevant columns and group by enhancer
comb_enh <- comb %>%
  select(grna_target, target_enh, is_target_enh, target_prom, is_target_prom, ctrl_perturbation, ctrl_type, enhancer_like_interaction, distToTSS_bin, trait) %>%
  distinct() %>%
  filter(ctrl_type != "promoter" | is.na(ctrl_type)) %>% #remove promoter elements
  filter(enhancer_like_interaction == TRUE) %>%
  filter(distToTSS_bin != "> 1Mb") %>%
  separate_rows(grna_target, sep = ",") %>%
  mutate(trait = case_when(trait == "Primary_biliary_cirrhosis" ~ "Primary_biliary_cholangitis",
                           TRUE ~ trait)) %>% #cirrhosis older name for same disease, rename and combine
  group_by(grna_target) %>%
  add_tally()

#List of enhancers to use for filtering of downstream results
# genes_by_enh <- comb_enh %>%
#   filter(!is.na(grna_target)) %>%
#   mutate(grna_target = str_replace_all(grna_target, "[:\\-]", "_")) %>%
#   group_by(grna_target) %>%
#   summarise(genes = list(target_prom), .groups = "drop") %>%
#   mutate(n_genes = lengths(genes)) %>%
#   filter(n_genes > 10) %>% #remove diseases where the gene list only contains 10 genes or less
#   select(-n_genes) %>%
#   deframe()

#Adjust list of enhancers and link to downstream genes
genes_by_enh_df <- comb_enh %>%
  filter(!is.na(grna_target)) %>%
  mutate(grna_target = str_replace_all(grna_target, "[:\\-]", "_")) %>%
  group_by(grna_target) %>%
  summarise(genes = list(sort(unique(target_prom))),.groups = "drop") %>%
  mutate(n_genes = lengths(genes)) %>%
  filter(n_genes > 10) %>%
  select(-n_genes)

#Group together enhancers that have exactly the same set of downstream genes
genes_by_enh_grouped <- genes_by_enh_df %>%
  mutate(gene_signature = map_chr(genes, ~ paste(.x, collapse = "|"))) %>%
  group_by(gene_signature) %>%
  summarise(
    enhancers = list(grna_target),
    genes = list(genes[[1]]),
    .groups = "drop"
  )

#List of enhancers to use for enrichment
genes_by_enh <- setNames(genes_by_enh_grouped$genes,
                         map_chr(genes_by_enh_grouped$enhancers, ~ paste(.x, collapse = "__")))

#genes_with_multiple_enh <- genes_by_enh[grepl("__", names(genes_by_enh2), fixed = TRUE)]

#Run enrichments
run_enrichment_sets(
  genes_by_enh,
  prefix = paste0(outdir, "enhancers"),
  gtf,
  all_genes_entrez,
  all_genes_ensembl,
  disease2gene,
  disease2name
)

#### SET OF GENES + DOWNSTREAM GENES ("HOP 2 GENES") LINKED TO ONE ENHANCER #### 

#select promoter screen results that are significant and with abs(logFC) > 0.2
res.filt <- res %>%
  filter(significant == TRUE) %>%
  filter(abs(logFC) > 0.2)
  
#Adjust list of enhancers and link to downstream genes (hop 1) + downstream genes (hop 2)
genes_by_enh_hop2_df <- comb_enh %>%
  filter(!is.na(grna_target)) %>%
  inner_join(res.filt, by = c("target_prom" = "perturbation"), relationship = "many-to-many") %>%
  select(grna_target, target_enh, is_target_enh, target_prom, is_target_prom, target_prom_hop2 = target, is_target_prom_hop2 = is_target,
         ctrl_perturbation, ctrl_type, enhancer_like_interaction, distToTSS_bin) %>%
  mutate(grna_target = str_replace_all(grna_target, "[:\\-]", "_")) %>%
  group_by(grna_target) %>%
  summarise(genes = list(sort(unique(c(target_prom), target_prom_hop2))),
            .groups = "drop") %>%
  mutate(n_genes = lengths(genes)) %>%
  filter(n_genes > 10) %>%
  select(-n_genes)

#Group together enhancers that have exactly the same set of downstream genes
genes_by_enh_hop2_grouped <- genes_by_enh_hop2_df %>%
  mutate(gene_signature = map_chr(genes, ~ paste(.x, collapse = "|"))) %>%
  group_by(gene_signature) %>%
  summarise(
    enhancers = list(grna_target),
    genes = list(genes[[1]]),
    .groups = "drop"
  )

#List of enhancers to use for enrichment
genes_by_enh_hop2 <- setNames(genes_by_enh_hop2_grouped$genes,
                         map_chr(genes_by_enh_hop2_grouped$enhancers, ~ paste(.x, collapse = "__")))


#Run enrichments
run_enrichment_sets(
  genes_by_enh_hop2,
  prefix = paste0(outdir, "enhancers_hop2"),
  gtf,
  all_genes_entrez,
  all_genes_ensembl,
  disease2gene,
  disease2name
)


#### SET OF PERTURBATIONS THAT CO-CLUSTER IN MOFA FACTORS ####

#MOFA perturbation clusters
genes_by_mofa_pert <- mofa_pert %>%
  mutate(leiden = paste0("MOFA_p_", leiden)) %>%
  group_by(leiden) %>%
  summarise(genes = list(contrast), .groups = "drop") %>%
  deframe()

#Run enrichments
run_enrichment_sets(
  genes_by_mofa_pert,
  prefix = paste0(outdir, "MOFA"),
  gtf,
  all_genes_entrez,
  all_genes_ensembl,
  disease2gene,
  disease2name
)

#MOFA downstream genes clusters
genes_by_mofa_ds <- mofa_ds %>%
  mutate(leiden = paste0("MOFA_ds_", leiden)) %>%
  group_by(leiden) %>%
  summarise(genes = list(variable), .groups = "drop") %>%
  deframe()

#Run enrichments
run_enrichment_sets(
  genes_by_mofa_ds,
  prefix = paste0(outdir, "MOFA"),
  gtf,
  all_genes_entrez,
  all_genes_ensembl,
  disease2gene,
  disease2name
)

