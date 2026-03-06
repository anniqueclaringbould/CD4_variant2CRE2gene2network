###############################################
## Script to perform differential expression ##
#### Using limma-trend on individual cells ####
############ Split by perturbation ############
###############################################

#### SETUP ####
#libraries
library(dplyr, quietly = T)
library(dtplyr, quietly = T)
library(tidyverse, quietly = T)
library(limma, quietly = T)
library(edgeR, quietly = T)
library(ggpubr, quietly = T)
library(Seurat, quietly = T)
library(BPCells, quietly = T)
library(Matrix, quietly = T)
library(data.table, quietly = T)

#directories
root <- "/g/steinmetz/project/otar_2063/promoter_screen/"
indir <- paste0(root, "interpretation/seurat_qc/output/combined/")
degdir <- paste0(root, "interpretation/DE/")
outdir <- paste0(degdir, "output/per_gene/")

#### READ DATA ####
seur <- readRDS(paste0(indir, "Seurat_merged_QC_2025-08-04.rds"))
guide_assignment <- readMM(paste0(root, "interpretation/guide_assignments/matrix/grna_assignment_matrix.mtx"))
guide_rows <- fread(paste0(root, "interpretation/guide_assignments/matrix/grna_assignment_matrix_rownames.txt"), header = F)
guide_cols <- fread(paste0(root, "interpretation/guide_assignments/matrix/grna_assignment_matrix_colnames.txt"), header = F)
rownames(guide_assignment) <- guide_rows$V1
colnames(guide_assignment) <- guide_cols$V1
nt_guides <- fread("/g/stegle/schrod/code/TCell/NT_guide_DEG_counts.csv") #list of non-targeting guides with the number of hits
covariates <- fread(paste0(root, "interpretation/cell_scoring/output/Representative_AUC_scores_per_cell.tsv")) #AUC scores used as covariates
genes_to_test <- fread("/g/stegle/schrod/code/TCell/gene_set.csv") #set of genes to test based on being highly variable + our targets

#### PREPARE ####

#add covariates to seurat meta data
seur <- AddMetaData(seur, covariates)

#prepare gRNA info

#convert to a long table of (cell, guide) for entries == 1
ga_dt <- data.table(
  guide = rownames(guide_assignment)[guide_assignment@i + 1],
  cell  = colnames(guide_assignment)[guide_assignment@j + 1]
)

#count guides assigned per cell and save perturbations
ga_dt <- lazy_dt(ga_dt) %>% #dtplyr pipeline (readable but fast)
  group_by(cell) %>%
  mutate(perturbation = sub("-[0-9]+$", "", guide),
         nr_gRNAs_assigned = n())

#save NT guides that have >15 differential expression effects
remove_guides <- nt_guides %>%
  filter(num_degs >= 15) %>%
  pull(V1)

#keep all cells with exactly one non-targeting guide since they should be included in each subset
#remove guides with too many off-target effects
ga_NT <- ga_dt %>%
  filter(nr_gRNAs_assigned == 1) %>%
  filter(perturbation == "GUIDE-NO-TARGET") %>%
  filter(!guide %in% remove_guides) %>%
  select(cell)

#keep all cells with exactly one targeting guide
ga_regular <- ga_dt %>%
  filter(nr_gRNAs_assigned == 1) %>%
  filter(perturbation != "GUIDE-NO-TARGET") %>%
  distinct(cell, perturbation)

#randomly split perturbations into groups
set.seed(1)
perturbation_groups <- ga_regular %>%
  ungroup() %>%
  select(perturbation) %>%
  distinct(perturbation) %>%
  mutate(group = sample(1:5, n(), replace = T))

# materialize once, before list + loop work
ga_dt <- as.data.table(ga_dt)
ga_NT <- as.data.table(ga_NT)
ga_regular <- as.data.table(ga_regular)
perturbation_groups <- as.data.table(perturbation_groups)

#make separate lists of cells
cell_ids <- perturbation_groups %>%
  left_join(ga_regular, by = "perturbation") %>%
  group_by(group) %>%
  group_split()  %>%
  lapply(function(df) {
    bind_rows(df, ga_NT) %>%
      pull(cell)
  })

#write list with cells for each split
cell_split_table <- tibble(group = rep(seq_along(cell_ids), lengths(cell_ids)),
                           cell  = unlist(cell_ids))

write.table(cell_split_table, paste0(outdir, "cell_splits.tsv"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# save count matrix (genes × cells)
message(format(Sys.time()), ': saving count matrix')
counts_mat <- seur@assays$RNA$counts
message(paste0(format(Sys.time()), ': Unfiltered count matrix has ', dim(counts_mat)[1], ' genes and ', dim(counts_mat)[2], ' cells'))

# save genes present in counts matrix
genes <- rownames(counts_mat)

#combine gene lists and filter count matrix
genes.keep <- unique(intersect(genes, genes_to_test$bpcells_name))
counts_mat <- counts_mat[genes.keep, , drop = FALSE]
message(paste0(format(Sys.time()), ': Filtered count matrix has ', dim(counts_mat)[1], ' genes and ', dim(counts_mat)[2], ' cells'))

#write matrix used for this analysis
write_rds(counts_mat, paste0(outdir, "Counts_matrix_limma_splits_", Sys.Date(), ".rds"))

#RUN LIMMA FOR EACH GROUP SEPARATELY
for (split in seq_along(cell_ids)) {
  #for (split in 1:2) {
  
  message(format(Sys.time()), ": Now running limma for split ", split)
  
  #define cells once, in the desired order
  cells_split <- cell_ids[[split]]
  
  #subset original table to those cells
  ga_split_dt <- ga_dt[cell %in% cells_split]
  
  #keep only cells with gRNA assignment in seurat object
  seur_sub <- subset(seur, cells = cells_split)
  
  #prepare metadata
  message(format(Sys.time()), ': Preparing meta data')
  
  #keep only cells with gRNA assignment in meta data
  meta <- seur_sub@meta.data %>%
    left_join(ga_split_dt, by = "cell") %>% 
    separate_wider_delim(guide, delim = "-", too_many = "merge", names = c("tmp", "tail")) %>% #split guide into prefix and the rest
    mutate(guide_nr = sub(".*-(\\d+)$", "\\1", tail),   #split tail into gene and guide_nr by LAST hyphen
           gene = sub("(.*)-\\d+$", "\\1", tail),
           guide = paste0(tmp, "-", gene)) %>%
    select(cell, gene, guide, donor, lane = sample_lane, sequencing_run = orig.ident, 
           percent.mt, percent.rb, nr_genes = nFeature_RNA, nr_UMIs = nCount_RNA,
           CellCycle.G2M, early, ICOS.CD38, Mito, CTLA4.CD38, Poor.Quality, iTreg_5d4h, Th1_5d) %>%
    slice(match(cells_split, cell)) %>% #restore cell order
    column_to_rownames("cell")
  
  # relevel so that non targeting "gene" is the baseline
  meta$gene <- relevel(factor(meta$gene), ref = "NO-TARGET")
  
  #keep only cells with gRNA assignment in counts matrix
  counts_mat_sub <- counts_mat[, cells_split, drop = FALSE]
  
  message(format(Sys.time()), ': Checking if cells are shared between count matrix and meta data')
  
  #check if cell order in meta and counts matrix are the same
  stopifnot(identical(colnames(counts_mat_sub), rownames(meta)))
  
  
  #### DEG DESIGN ####
  
  message(format(Sys.time()), ': Starting design matrix')
  
  #design matrix
  design <- sparse.model.matrix(~ gene + sequencing_run + 
                                  nr_UMIs + percent.mt + percent.rb + 
                                  CellCycle.G2M + early + ICOS.CD38 + 
                                  Mito + CTLA4.CD38 + Poor.Quality + iTreg_5d4h, 
                                data = meta) 
  
  message(format(Sys.time()), ': Design matrix dimensions are: ', dim(design)[1], ' by ', dim(design)[2])
  
  #create DGEList object
  suppressWarnings(dge <- DGEList(counts = counts_mat_sub))
  
  #calculate size factors
  message(format(Sys.time()), ': Calculating norm factors')
  dge <- calcNormFactors(dge)
  
  #normalize to counts per million
  message(format(Sys.time()), ': Normalizing to CPM')
  logCPM <- edgeR::cpm(dge, log=TRUE, prior.count = 3)
  
  #run differential expression for all groups (genes)
  message(format(Sys.time()), ': Running DE fit')
  fit <- lmFit(logCPM, design)
  message(format(Sys.time()), ': Running eBayes fit')
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  #write output for each gene separately
  coef_names <- colnames(fit$coefficients) # get all coefficient names
  gene_coefs <- coef_names[grepl("^gene", coef_names)] # keep only those starting with "gene"
  
  # loop through each gene for writing the output
  message(format(Sys.time()), ': Printing per gene output')
  
  for (coef in gene_coefs) {
    tt <- topTreat(fit, coef = coef, n = Inf, adjust.method = "BH")
    outfile <- file.path(outdir, paste0("DEG_single_cell_limma_ebayes_split_", split, "_NT_removed_all_covariates_prefiltered_genes_", coef, "_", Sys.Date(), ".tsv"))
    fwrite(tt, outfile, sep = "\t", row.names = TRUE) 
  }
  
  message(format(Sys.time()), ": Finished calculating differential expression for split ", split)
  
}

