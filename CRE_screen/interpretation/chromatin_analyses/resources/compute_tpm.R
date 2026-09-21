## Compute TPM values from gene expression matrix. This is technically computing-counts-per-million
## (CPM) based on 3' UMI counts. Since each UMI represents one molecule (i.e. transcript), this
## calculation is equivalent to transcripts-per-million (TPM) 

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(rtracklayer)
  library(Matrix)
})

# required input files
annot_file <- "/g/steinmetz/project/otar/forAndreas/hg38.gtf"
dge_matrix_file <- "/g/steinmetz/project/otar_2063/enhancer_screen/interpretation/chromatin_analyses/resources/dge.txt.gz"

# load genome annotations and extract unique gene ids and gene symbols
annot <- import(annot_file, format = "gtf")
genes <- annot[annot$type == "gene"]
gene_ids <- mcols(genes) %>% 
  as.data.frame() %>% 
  select("gene_id", "gene_name") %>% 
  mutate(gene_id = sub("\\..+", "", gene_id)) %>% 
  distinct()

# load file containing 10x UMI matrix and create sparse matrix
dge <- fread(dge_matrix_file)
gene_symbols <- dge$GENE
dge <- as(data.matrix(dge[, -1]), "sparseMatrix")
rownames(dge) <- gene_symbols

# compute TPM for each gene
txs_per_gene <- rowSums(dge)
tpm <- txs_per_gene / sum(txs_per_gene) * 1e6

# create data frame containing TPM and add gene ids
tpm <- data.table(gene = names(tpm), tpm = tpm)
tpm <- left_join(tpm, gene_ids, by = c("gene" = "gene_name"))

# save to output file
fwrite(tpm, file = "gene_tpm.csv.gz")
