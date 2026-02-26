#%%
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np

import os
import sys
import scipy.io
import click

import seaborn as sns
import seaborn.objects as so
import matplotlib.pyplot as plt
from plotnine import *

import logging
import click
import h5py
from anndata._io.specs import read_elem


logging.basicConfig(format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s", stream=sys.stdout)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

def get_diff_expression_results(adata):
    df_list = []
    for var in ["names", "scores", "pvals", "pvals_adj", "logfoldchanges"]:
        tmp_df = pd.DataFrame(adata.uns["rank_genes_groups"][var])
        tmp_df = tmp_df.melt(var_name="group", value_name=var)
        df_list.append(tmp_df)
    df = pd.concat(df_list, axis=1)
    df = df.loc[:, ~df.columns.duplicated()]
    return df


"""
start_idx = 0
end_idx = 5
target_column = "target"
control_group = "NO-TARGET"
h5ad_file = f"/g/stegle/schrod/data/T_Cell/promoter_single_guide_hvg.h5ad"
results_dir = "/g/stegle/schrod/code/TCell/results/DEG/"
"""



@click.command()
@click.option('-s', '--start-idx', type=int, required=True, help='Start of perturbation range')
@click.option('-e', '--end-idx', type=int, required=True, help='End of perturbation range')
@click.option('-tc', '--target-column', default="target", help='Column with target names')
@click.option('-tcn', '--target-column-idx', default="target_idx", help='Column with target idxs')
@click.option('-c', '--control-group', default="NO-TARGET", help='Control group identifier')
@click.option('-adata-path', '--adata-path', default="/g/stegle/schrod/data/T_Cell/promoter_single_guide_hvg.h5ad", help='Adata path')
@click.option('-results-dir', '--results-dir', default="/g/stegle/schrod/code/TCell/results/DEG/", help='Result directory')
@click.option('-cutoff', '--cutoff', default=5, help='min number of replicates/cells')
@click.option('-nt_cutoff', '--nt_cutoff', type=int, default=-1, help='cutoff of num DGEs per NT guide')
@click.option('-donor', '--donor', type=str, default="all", help='To subset by donor')
@click.option('-remove_donor', '--remove_donor', type=str, default="keep_all", help='To subset by donor')
@click.option('-num_guides', '--num_guides', type=int, default=0, help='Set number of NT guides')
@click.option('-deg_test', '--deg_test', type=str, default="wilcoxon", help='Set test')
def compute_partial_DEGs(start_idx, end_idx, target_column, target_column_idx, control_group, adata_path, results_dir, cutoff, nt_cutoff, donor, remove_donor, num_guides, deg_test):
    """
    Calculate differential expressed genes between numerical start and end indices to facilitate parallel computation
    """
    logger.info(f"Read adata from: {adata_path}")
    adata = ad.read_h5ad(os.path.join(adata_path))

    if donor != "all":
        logger.info(f": Subset to donor {donor}")
        adata = adata[adata.obs.donor == donor]
    
    if remove_donor != "keep_all":
        logger.info(f": Remove donor {remove_donor}")
        adata = adata[adata.obs.donor != remove_donor]


    if nt_cutoff > 0:
        logger.info(f": Remove NT guides with less than {nt_cutoff} DEGs")
        nt_guide_counts = pd.read_csv("/g/stegle/schrod/code/TCell/NT_guide_DEG_counts.csv", index_col=0)
        all_nt_guides = list(adata.obs.guide[adata.obs.target=="NO-TARGET"].unique())
        nt_guides_to_keep = nt_guide_counts[nt_guide_counts.num_degs < nt_cutoff].index.tolist()

        if num_guides > 0:
            nt_guides_to_keep = nt_guides_to_keep[:num_guides]

        remove_nt_guides = list(set(all_nt_guides) - set(nt_guides_to_keep))
        adata = adata[~adata.obs["guide"].isin(remove_nt_guides)]
    
    if target_column == "guide":
        tmp = adata.obs[target_column].astype(str)
        tmp[adata.obs["target"] == "NO-TARGET"] = "NO-TARGET"
        adata.obs[target_column] = tmp.astype("category")
        logger.info(f": Converted guide column to categorical with NO-TARGET as control group")

    logger.info(f": Remove targets with less than {cutoff} cells")
    adata = adata[adata.obs.groupby(target_column)[target_column].transform('size')>cutoff]

    idxs = np.arange(start_idx, end_idx)
    mask = adata.obs[target_column_idx].isin(idxs)
    target_list = list(adata.obs[target_column][mask].unique())
    
    all_results = []
    for target in target_list:
        logger.info(f"Processing target: {target}")

        # Determine nonzero genes in pseudobulk for this target
        adata_sub = adata[adata.obs[target_column] == target]
        gene_sums = np.array(adata_sub.X.sum(0)).ravel()
        nonzero_genes = adata_sub.var_names[gene_sums != 0]

        if len(nonzero_genes) == 0:
            logger.info(f"Skipping target {target}: all genes are 0.")
            continue

        if target == control_group:
            logger.info(f"Skipping target {target}: is control group.")
            continue

        # Subset experimental AnnData to control + target cells, keeping only nonzero genes
        adata_sub = adata[adata.obs[target_column].isin([target, control_group])]
        adata_sub = adata_sub[:, adata_sub.var_names.isin(nonzero_genes)].copy()

        # Run DE only on nonzero genes
        sc.tl.rank_genes_groups(
            adata_sub,
            groupby=target_column,
            groups=[target],
            reference=control_group,
            method=deg_test
        )

        df_res = get_diff_expression_results(adata_sub)
        all_results.append(df_res)

    if len(all_results) == 0:
        logger.info("No DEGs found for selected targets.")
        return

    df_all = pd.concat(all_results, ignore_index=True)

    os.makedirs(results_dir, exist_ok=True)
    out_path = os.path.join(results_dir, f"DEG_filtered_from_{start_idx}_to_{end_idx-1}.csv")
    df_all.to_csv(out_path, index=False)

    logger.info(f"Saved filtered DEGs to {out_path}")

if __name__ == "__main__":
    compute_partial_DEGs()
