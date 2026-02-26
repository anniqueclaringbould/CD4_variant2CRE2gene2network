# %%
import pandas as pd
import anndata as ad
import scanpy as sc
import numpy as np

import os
import sys
import click

os.chdir("/g/stegle/schrod/code/TCell")
from utils import *

logger = get_logger(__name__)


def get_diff_expression_results(adata):
    df_list = []
    for var in ["names", "scores", "pvals", "pvals_adj", "logfoldchanges"]:
        tmp_df = pd.DataFrame(adata.uns["rank_genes_groups"][var])
        tmp_df = tmp_df.melt(var_name="group", value_name=var)
        df_list.append(tmp_df)
    df = pd.concat(df_list, axis=1)
    df = df.loc[:, ~df.columns.duplicated()]
    return df


@click.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.option('--iteration', type=int, required=True, help='Iteration number (0-based)')
@click.option('--n-guides', type=int, default=3, help='Number of guides to select for Group A')
@click.option('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/wilcoxon_calibration_NT", help='Output directory for results')
@click.option('--h5ad-file', type=str, default="/g/stegle/schrod/data/T_Cell/promoter_bulk_profiles_1guide_non_targeting.h5ad", help='Path to input h5ad file')
@click.option('--min-cells', type=int, default=5, help='Minimum cells per target group')
def main(iteration, n_guides, output_dir, h5ad_file, min_cells):

    logger.info(f"=== Iteration {iteration} ===")
    logger.info(f"Parameters: n_guides={n_guides}, min_cells={min_cells}")

    # Load data
    logger.info(f"Loading anndata object from: {h5ad_file}")
    adata = ad.read_h5ad(h5ad_file)
    adata.obs["target"] = adata.obs.guide.apply(lambda x: "-".join(x.split("-")[1:-1]))
    logger.info(f"Found {adata.obs.target.nunique()} unique targets.")

    # Get NT cells only
    adata_nt = adata[adata.obs['target'] == 'NO-TARGET'].copy()
    logger.info(f"Number of NO-TARGET cells: {adata_nt.n_obs}")

    # Get unique NT guides
    nt_guides = adata_nt.obs.guide.unique()
    n_total_guides = len(nt_guides)
    logger.info(f"Total NT guides: {n_total_guides}")

    if n_guides >= n_total_guides:
        logger.error(f"Cannot select {n_guides} guides from {n_total_guides} total guides")
        sys.exit(1)

    # Randomly select guides for Group A
    np.random.seed(iteration)
    shuffled_guides = np.random.permutation(nt_guides)

    group_A_guides = shuffled_guides[:n_guides]
    group_B_guides = shuffled_guides[n_guides:]

    logger.info(f"Group A: {len(group_A_guides)} guides, Group B: {len(group_B_guides)} guides")
    logger.info(f"Group A guides: {', '.join(group_A_guides)}")

    # Assign groups
    adata_nt.obs['target'] = adata_nt.obs.guide.apply(
        lambda x: 'targeting' if x in group_A_guides else 'NO-TARGET'
    )

    n_targeting = (adata_nt.obs.target == 'targeting').sum()
    n_control = (adata_nt.obs.target == 'NO-TARGET').sum()
    logger.info(f"Group A cells: {n_targeting}")
    logger.info(f"Group B cells: {n_control}")

    # Filter groups with too few cells
    if n_targeting < min_cells:
        logger.error(f"Group A has only {n_targeting} cells (< {min_cells})")
        sys.exit(1)
    if n_control < min_cells:
        logger.error(f"Group B has only {n_control} cells (< {min_cells})")
        sys.exit(1)

    # Determine nonzero genes in targeting group
    adata_targeting = adata_nt[adata_nt.obs.target == 'targeting']
    gene_sums = np.array(adata_targeting.X.sum(0)).ravel()
    nonzero_genes = adata_targeting.var_names[gene_sums != 0]
    logger.info(f"Nonzero genes in Group A: {len(nonzero_genes)}")

    # Subset to nonzero genes
    adata_nt = adata_nt[:, adata_nt.var_names.isin(nonzero_genes)].copy()

    # Run Wilcoxon rank-sum test
    logger.info("Running Wilcoxon rank-sum test: targeting vs NO-TARGET...")
    sc.tl.rank_genes_groups(
        adata_nt,
        groupby='target',
        groups=['targeting'],
        reference='NO-TARGET',
        method='wilcoxon'
    )

    result_df = get_diff_expression_results(adata_nt)

    # Add iteration info
    result_df['iteration'] = iteration
    result_df['n_guides_selected'] = n_guides
    result_df['group_A_guides'] = ','.join(group_A_guides)

    # Save results
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, f"iteration_{iteration:03d}.csv")
    result_df.to_csv(output_path, index=False)
    logger.info(f"Results saved to: {output_path}")

    # Log summary statistics
    n_sig_005 = (result_df['pvals_adj'] < 0.05).sum()
    n_sig_01 = (result_df['pvals_adj'] < 0.1).sum()
    logger.info(f"Significant genes (p<0.05): {n_sig_005}")
    logger.info(f"Significant genes (p<0.1): {n_sig_01}")

    logger.info(f"Iteration {iteration} complete!")

if __name__ == "__main__":
    main()
