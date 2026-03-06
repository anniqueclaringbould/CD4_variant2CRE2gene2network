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
outdir <- paste0(degdir, "output/per_gene/non_targeting/")

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

#keep all cells with exactly one non-targeting guide since they should be the main dataset here
#remove guides with too many off-target effects
ga_NT <- ga_dt %>%
  filter(nr_gRNAs_assigned == 1) %>%
  filter(perturbation == "GUIDE-NO-TARGET") %>%
  select(cell, guide)

#randomly assign NT guides into groups of perturbations
set.seed(1)

perturbation_groups <- ga_NT %>%
  ungroup() %>%
  distinct(guide) %>%
  slice_sample(prop = 1) %>% # shuffle guides
  mutate(group = rep(seq_len(ceiling(n() / 3)), each = 3)[seq_len(n())]) #assign each 3 guides to a group

# materialize once, before list + loop work
ga_dt <- as.data.table(ga_dt)
ga_NT <- as.data.table(ga_NT)
perturbation_groups <- as.data.table(perturbation_groups)

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
write_rds(counts_mat, paste0(outdir, "Counts_matrix_limma_fake_perturbations_", Sys.Date(), ".rds"))

#RUN LIMMA FOR EACH FAKE PERTURBATION (3 GUIDES) SEPARATELY
groups <- sort(unique(perturbation_groups$group))

for (g in groups) {
  
  message(format(Sys.time()), ": Running limma for fake perturbation group ", g)
  
  ## cells in this fake perturbation
  cells_group <- perturbation_groups %>%
    filter(group == g) %>%
    left_join(ga_NT, by = "guide") %>%
    pull(cell) %>%
    unique()
  
  ## all other NT cells
  cells_other <- perturbation_groups %>%
    filter(group != g) %>%
    left_join(ga_NT, by = "guide") %>%
    pull(cell) %>%
    unique()
  cells_use <- c(cells_group, cells_other)
  
  ## subset Seurat and reorder counts
  seur_sub <- subset(seur, cells = cells_use)
  counts_mat_sub <- counts_mat[, cells_use, drop = FALSE]
  
  ## prepare metadata
  meta <- seur_sub@meta.data %>%
    mutate(
      cell = rownames(.),
      test_group = ifelse(cell %in% cells_group, "group", "other"),
      test_group = relevel(factor(test_group), ref = "other")
    ) %>%
    select(cell, test_group, donor, lane = sample_lane, sequencing_run = orig.ident,
           percent.mt, percent.rb, nr_genes = nFeature_RNA, nr_UMIs = nCount_RNA,
           CellCycle.G2M, early, ICOS.CD38, Mito, CTLA4.CD38, Poor.Quality, iTreg_5d4h) %>%
    slice(match(cells_use, cell)) 
  
  stopifnot(identical(colnames(counts_mat_sub), rownames(meta)))
  
  ## design matrix
  design <- sparse.model.matrix(
    ~ test_group + sequencing_run +
      nr_UMIs + percent.mt + percent.rb +
      CellCycle.G2M + early + ICOS.CD38 +
      Mito + CTLA4.CD38 + Poor.Quality + iTreg_5d4h,
    data = meta
  )
  
  ## limma
  dge <- DGEList(counts = counts_mat_sub)
  dge <- calcNormFactors(dge)
  logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 3)
  
  fit <- lmFit(logCPM, design)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  
  ## extract result
  tt <- topTreat(fit, coef = "test_groupgroup", n = Inf, adjust.method = "BH")
  
  outfile <- file.path(
    outdir,
    paste0("DEG_single_cell_limma_fakePerturbation_", g,
           "_vs_all_", Sys.Date(), ".tsv")
  )
  fwrite(tt, outfile, sep = "\t", row.names = TRUE)
  
  message(format(Sys.time()), ": Finished fake perturbation group ", g)
}






