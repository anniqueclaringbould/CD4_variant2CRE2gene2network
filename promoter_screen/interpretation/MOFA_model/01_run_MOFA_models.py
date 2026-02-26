import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
import muon as mu
import mofaflex as mfl

def load_deg_results(file_path, pval_threshold = None):
    df_list = []
    for file_name in os.listdir(file_path):
        if file_name == "DEG_full.csv":
            continue
        df_tmp = pd.read_csv(os.path.join(file_path, file_name), index_col = 0)
        if pval_threshold is not None:
            df_tmp = df_tmp[df_tmp["adj_p_value"] < pval_threshold]
        df_list.append(df_tmp)
    df1 = pd.concat(df_list)
    return df1

# read in data (all LFCs instead of only significant ones)
lfcs = load_deg_results('data/pseudobulk_deseq2_chunks50_min3bulks_ribosomal/', pval_threshold = 0.2)

# remove on-target effect such that they don't distort clustering
lfcs = lfcs[lfcs['variable'] != lfcs['contrast']]

# filter for gRNA targets and response genes that are at least significant 10 times
lfcs_sig = load_deg_results('data/pseudobulk_deseq2_chunks50_min3bulks_ribosomal/', pval_threshold = 0.1)
grna_targets = lfcs_sig['contrast'].unique()
counts = lfcs_sig['contrast'].value_counts()
names = counts[counts > 10].index.tolist()
print('Number of target genes with at least one significant DE gene:', len(grna_targets))
print('Number of target genes with more than 10 significant DE genes:', len(names))
lfcs = lfcs[lfcs['contrast'].isin(names)]

response_genes = lfcs_sig['variable'].unique()
counts = lfcs_sig['variable'].value_counts()
response_names = counts[counts > 10].index.tolist()
print('Number of response genes that are significant in at least one perturbation:', len(response_genes))
print('Number of response genes that are significant in more than 10 selected perturbations:', len(response_names))
lfcs = lfcs[lfcs['variable'].isin(response_names)]

lfcs_wide = lfcs.pivot(index="contrast", columns="variable", values="log_fc")
adata = ad.AnnData(lfcs_wide)
mdata = mu.MuData({"gene_expression": adata})

# run MOFA for various seeds with various number of factors
for n_factors in [5, 10, 15, 20, 25, 30]:
    for run in range(10):
        mofa_model = mfl.MOFAFLEX(
            mdata,
            mfl.ModelOptions(n_factors=n_factors, weight_prior = "Horseshoe", factor_prior = "Horseshoe", likelihoods="Normal"),
            mfl.TrainingOptions(device = "cuda", seed=42 + run, 
                                save_path='models/mofa_lr_n' + str(n_factors) + "_" + str(run), lr=0.005, 
                                early_stopper_patience=1000, max_epochs = 20000), 
        )

