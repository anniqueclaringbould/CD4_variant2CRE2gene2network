import numpy as np
import scanpy as sc
import pandas as pd
import anndata as ad
import logging
import matplotlib.pyplot as plt
import sys
import os
from tqdm import tqdm

def get_logger(name=__name__):
    logging.basicConfig(format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s", stream=sys.stdout)
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    return logger

def norm_log(adata, logger=None):
    if logger:
        logger.info(f"Normalize data...")
    sc.pp.normalize_total(adata)

    if logger:
        logger.info(f"Log transform data...")
    sc.pp.log1p(adata)

def preprocess(adata, logger=None):
    if logger:
        logger.info(f"Compute pca...")
    sc.tl.pca(adata)
    
    if logger:
        logger.info(f"Compute nearest neighbours...")
    sc.pp.neighbors(adata)

    if logger:
        logger.info(f"Compute umap encoding...")
    sc.tl.umap(adata)

def plot_umap(adata, key, figure_path):
    fig, ax = plt.subplots() 
    sc.pl.umap(adata, color=key, ax=ax)
    fig.tight_layout()
    fig.savefig(figure_path, dpi=500)

def get_diff_expression_results(adata):
    df_list = []
    # loop over entries of the adata.uns["rank_gens_groups"] slot and write to dataframe
    for var in ["names", "scores", "pvals", "pvals_adj", "logfoldchanges"]:
        tmp_df = pd.DataFrame(adata.uns["rank_genes_groups"][var])
        # add one column indicating the comparison group 
        tmp_df = tmp_df.melt(var_name="group", value_name=var)
        df_list.append(tmp_df)
    df = pd.concat(df_list, axis=1)
    df = df.loc[:, ~df.columns.duplicated()]
    return df

def load_deg_results(file_path, pval_threshold = None):
    df_list = []
    for file_name in tqdm(os.listdir(file_path)):
        if file_name == "DEG_full.csv":
            continue
        df_tmp = pd.read_csv(os.path.join(file_path, file_name))
        if pval_threshold is not None:
            df_tmp = df_tmp[df_tmp["pvals_adj"] < pval_threshold]
        df_list.append(df_tmp)
    df1 = pd.concat(df_list)
    return df1


def median_of_ratios(counts: np.ndarray):
    """
    DESeq2 Median of Ratios normalization (pure NumPy).
    counts: genes x samples
    """
    mask_nonzero = counts > 0
    geom_means = np.exp(
        np.sum(np.log(counts, where=mask_nonzero), axis=1) /
        np.maximum(mask_nonzero.sum(axis=1), 1)
    )
    geom_means[geom_means == 0] = np.nan
    ratios = counts / geom_means[:, np.newaxis]
    size_factors = np.nanmedian(ratios, axis=0)
    norm_counts = counts / size_factors[np.newaxis, :]
    return size_factors, norm_counts


def get_efficiancy(adata, target_column = "perturbed_gene_name", no_target_key = "NTC"):
    all_targets = adata.var_names[adata.var_names.isin(adata.obs[target_column].unique())]
    results = {}
    adata_control = adata[adata.obs[target_column] == no_target_key]
    control_means = np.array(adata_control.X.mean(0)).flatten()

    for gene in tqdm(all_targets):
        adata_target = adata[adata.obs[target_column] == gene, gene]
        pert_mean = adata_target.X.mean()
        control_mean = control_means[np.argwhere(adata.var_names==gene)[0]]
        
        results[gene] = {
            'control_mean': control_mean,
            'perturbed_mean': pert_mean,
            'fraction': pert_mean / control_mean
        }

    df_results = pd.DataFrame(results).T
    df_results = df_results.reset_index()
    df_results.columns = ['gene', 'control_mean', 'perturbed_mean', 'fraction']

    df_results["fraction"] = df_results["fraction"].astype("float")
    df_results["control_mean"] = df_results["control_mean"].astype("float")

    
    return df_results