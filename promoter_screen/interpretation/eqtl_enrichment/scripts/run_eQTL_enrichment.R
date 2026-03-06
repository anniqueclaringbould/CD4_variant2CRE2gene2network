########################################
## Script to quantify eQTL enrichment ##
### cis-trans links from eQTLGen p2 ####
########################################

#### SETUP ####
#libraries
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(rtracklayer)

#directories
root    <- "/g/steinmetz/project/otar_2063/"
eqtldir <- paste0(root, "promoter_screen/interpretation/eqtl_enrichment/")
outdir  <- paste0(eqtldir, "output/")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)


#### INPUT — choose method here ####

method <- "pseudobulk_deseq"     # "pseudobulk_deseq" | "pseudobulk" | "sceptre" | "sceptre_empirical" | "limma" | "wilcox"

limma_file  <- paste0(root, "promoter_screen/interpretation/DE/output/combined/DEG_single_cell_limma_ebayes_split_NT_removed_all_covariates_prefiltered_genes_gene_combined_2026-02-03.tsv.gz")
#sceptre_file <- paste0(root, "promoter_screen/interpretation/diff_expr_analyses/process_trans_analysis/sceptre_pvalues/all_pairs.tsv.gz")
sceptre_file <- paste0(root, "promoter_screen/interpretation/diff_expr_analyses/process_trans_analysis/sceptre_pvalues_filt_genes/all_pairs.tsv.gz")
sceptre_emp_file <- paste0(root, "promoter_screen/interpretation/diff_expr_analyses/process_trans_analysis/gene_level_empirical_pvalues/all_pairs_empirical_pvalues_per_gene.tsv.gz")
pbulk_files <- list.files("/g/stegle/schrod/code/TCell/results_pseudobulk/DEG_bulk_rerun/", full.names = TRUE)
#pbulk_deseq_files <- list.files("/g/stegle/schrod/code/TCell/pseudobulk_deseq2_chunks50/", full.names = TRUE)
pbulk_deseq_files <- list.files("/g/stegle/schrod/code/TCell/pseudobulk_deseq2_chunks50_min3bulks_ribosomal/", full.names = TRUE)
wilcox_files <- list.files("/g/stegle/schrod/code/TCell/results_single_cell/DEG_15_nt_cutoff/", full.names = TRUE)

pval_threshold <- 0.1

#### LOADERS ####

load_results <- function(method){
  
  if(method == "limma"){
    dm <- fread(limma_file)
    
    list(
      res = dm %>%
        transmute(
          perturbation   = Perturbation,
          target = Target,
          significant = adj.P.Val < pval_threshold
        ),
      genes = unique(dm$Target)
    )
    
  } else if(method == "wilcox"){
    dm <- rbindlist(lapply(wilcox_files, fread))
    
    list(
      res = dm %>%
        transmute(
          perturbation   = group,
          target = names,
          significant = pvals_adj < pval_threshold
        ),
      genes = unique(dm$names)
    )
    
  } else if(method == "pseudobulk"){
    dm <- rbindlist(lapply(pbulk_files, fread))
    
    list(
      res = dm %>%
        transmute(
          perturbation   = group,
          target = names,
          significant = pvals_adj < pval_threshold
        ),
      genes = unique(dm$names)
    )
    
  } else if(method == "pseudobulk_deseq"){
    dm <- rbindlist(lapply(pbulk_deseq_files, fread))
    
    list(
      res = dm %>%
        transmute(
          perturbation   = contrast,
          target = variable,
          significant = adj_p_value < pval_threshold
        ),
      genes = unique(dm$contrast)
    )
    
  } else if(method == "sceptre"){
    dm <- fread(sceptre_file)
    
    list(
      res = dm %>%
        transmute(
          perturbation   = grna_target,
          target = gene_symbol,
          significant = p_value_adj < pval_threshold
        ),
      genes = unique(dm$gene_symbol)
    )
    
  } else if(method == "sceptre_empirical"){
    dm <- fread(sceptre_emp_file)
    
    list(
      res = dm %>%
        transmute(
          perturbation   = grna_target,
          target = gene_symbol,
          significant = pval_adj_per_gene < pval_threshold
        ),
      genes = unique(dm$gene_symbol)
    )
  }
}


#### READ SHARED DATA ####

eqtls <- fread(paste0(eqtldir, "input/CisTransRelationshipsBasedOnCololocization2025-11-14.txt"))

gtf <- readGFF(paste0(root, "references/alignment_reference/hg38/gencode_gtf/default/hg38.gtf")) %>%
  mutate(gene_id = gsub("\\..*|_.*", "", gene_id)) %>%
  select(gene_id, gene_name) %>%
  distinct()

#function to convert ENSEMBL IDs to gene IDs
ensembl_to_symbols <- function(x){
  tibble(gene_id = x) %>%
    left_join(gtf, by = "gene_id") %>%
    pull(gene_name)
}


#### RUN ANALYSIS ####

loaded <- load_results(method)

res   <- loaded$res
genes <- loaded$genes

eqtls.filt <- eqtls %>%
  mutate(
    CisGene   = ensembl_to_symbols(gsub("\\..*", "", CisGene)),
    TransGene = ensembl_to_symbols(gsub("\\..*", "", TransGene))
  ) %>%
  filter(CisGene %in% genes, TransGene %in% genes)

res_ann <- res %>%
  left_join(eqtls.filt, by = c("perturbation" = "CisGene",
                          "target" = "TransGene"))

cont.table <- table(
  eqtl_annotated = !is.na(res_ann$Direction),
  significant   = res_ann$significant
)

fwrite(
  as.data.frame(cont.table),
  paste0(outdir, "contingency_table_gene-gene_pairs_", method, ".tsv"),
  sep = "\t"
)

print(fisher.test(cont.table))
print(chisq.test(cont.table))

res_sig <- res_ann %>%
  filter(!is.na(Direction), significant)

fwrite(
  res_sig,
  paste0(outdir, "promoter_screen_with_overlapping_cis-trans-eQTLs_", method, ".tsv"),
  sep = "\t"
)
