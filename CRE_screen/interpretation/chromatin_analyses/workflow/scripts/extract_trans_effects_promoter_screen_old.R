## Extract trans-acting effects from promoter screen results

# save.image("RDA/extractTransProm.rda")
# stop()

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(readr)
})

# set seed for random sampling
set.seed(snakemake@params$seed)

# column schema for parquet files from trans analysis
sceptre_schema <- schema(
  response_id = string(),
  grna_target = string(),
  n_nonzero_trt = int32(),
  n_nonzero_cntrl = int32(),
  pass_qc = boolean(),
  p_value = double(),
  log_2_fold_change = double()
)

# directory containing results files
results_dir <- unique(dirname(snakemake@input$results))
if (length(results_dir) > 1) {
  stop("More than one directory for screen results provided", call. = FALSE)
}

# trans analysis results
ds <- open_dataset(results_dir, schema = sceptre_schema)

# load gene symbols for each response_id
genes <- read_tsv(snakemake@input$features, show_col_types = FALSE)

# get all trans perturbation-gene pairs based on perturbation and target gene names
trans_pairs <- ds %>% 
  left_join(., genes, by = "response_id") %>% 
  filter(grna_target != gene_symbol) %>% 
  compute()

# compute FDR on a sample of the data
fdr <- trans_pairs %>% 
  slice_sample(n = 5e6) %>% 
  collect() %>% 
  distinct() %>% 
  mutate(p_value_adj = p.adjust(p_value, method = "fdr"))

# get nominal p-value threshold corresponding to 10% FDR
pval_threshold <- fdr %>% 
  filter(p_value_adj < 0.1) %>% 
  pull(p_value) %>% 
  max()

# filter for all significant trans hits
sig_trans_pairs <- trans_pairs %>% 
  filter(p_value <= pval_threshold) %>% 
  collect() %>% 
  distinct()

# reformat for output
sig_trans_pairs <- sig_trans_pairs %>% 
  select(promoter = grna_target, trans_gene = gene_symbol) %>% 
  arrange(promoter, trans_gene)

# save to output file
write_csv(sig_trans_pairs, file = snakemake@output[[1]])
