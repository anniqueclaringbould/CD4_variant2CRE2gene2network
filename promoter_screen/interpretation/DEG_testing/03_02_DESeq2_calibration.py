# %%
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
import scipy.sparse as sp
import itertools
from tqdm import tqdm
from pertpy.tools._differential_gene_expression._pydeseq2 import PyDESeq2

os.chdir("/g/stegle/schrod/code/TCell")
from utils import *

logger = get_logger(__name__)

@click.command(context_settings=dict(ignore_unknown_options=True, allow_extra_args=True))
@click.option('--iteration', type=int, required=True, help='Iteration number (0-based)')
@click.option('--n-guides', type=int, default=3, help='Number of guides to select for Group A')
@click.option('--n-cpus-fit', type=int, default=20, help='Number of CPUs for model fitting')
@click.option('--n-cpus-test', type=int, default=20, help='Number of CPUs for contrast testing')
@click.option('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/pseudobulk_deseq2_NT", help='Output directory for results')
@click.option('--h5ad-file', type=str, default="/g/stegle/schrod/data/T_Cell/DESeq2_donor_guide_bulks.h5ad", help='Path to input h5ad file')
@click.option('--min-cells', type=int, default=5, help='Minimum cells per pseudobulk')
def main(iteration, n_guides, n_cpus_fit, n_cpus_test, output_dir, h5ad_file, min_cells):
    """
    Run DESeq2 on a random split of NT guides for a single iteration.
    
    Example usage:
        python 07_05_DESeq2_NTC_slurm.py --iteration 0
        python 07_05_DESeq2_NTC_slurm.py --iteration 1 --n-guides 5
    """
    
    logger.info(f"=== Iteration {iteration} ===")
    logger.info(f"Parameters: n_guides={n_guides}, min_cells={min_cells}")
    logger.info(f"CPUs: fit={n_cpus_fit}, test={n_cpus_test}")
    
    # Load data
    logger.info(f"Loading bulked anndata object from: {h5ad_file}")
    bulks = sc.read_h5ad(h5ad_file)
    bulks.obs["target"] = bulks.obs.guide.apply(lambda x: "-".join(x.split("-")[1:-1]))
    logger.info(f"Found {bulks.obs.target.nunique()} unique targets.")
    
    # Filter by cell count
    bulks = bulks[bulks.obs.n_cells >= min_cells]
    logger.info(f"Filtered to {bulks.n_obs} pseudobulks with >= {min_cells} cells")
    
    # Get NT bulks only
    bulks_nt = bulks[bulks.obs['target'] == 'NO-TARGET']
    logger.info(f"Number of NO-TARGET bulks: {bulks_nt.n_obs}")
    
    # Get unique NT guides
    nt_guides = bulks_nt.obs.guide.unique()
    n_total_guides = len(nt_guides)
    logger.info(f"Total NT guides: {n_total_guides}")
    
    if n_guides >= n_total_guides:
        logger.error(f"Cannot select {n_guides} guides from {n_total_guides} total guides")
        sys.exit(1)
    
    # Randomly select guides for Group A
    np.random.seed(iteration)  # For reproducibility
    shuffled_guides = np.random.permutation(nt_guides)
    
    group_A_guides = shuffled_guides[:n_guides]  # Select n_guides
    group_B_guides = shuffled_guides[n_guides:]  # All remaining guides
    
    logger.info(f"Group A: {len(group_A_guides)} guides, Group B: {len(group_B_guides)} guides")
    logger.info(f"Group A guides: {', '.join(group_A_guides)}")
    
    # Create a copy and assign groups
    bulks_split = bulks_nt.copy()
    bulks_split.obs['target'] = bulks_split.obs.guide.apply(
        lambda x: 'targeting' if x in group_A_guides else 'NO-TARGET'
    )
    
    logger.info(f"Group A bulks: {(bulks_split.obs.target == 'targeting').sum()}")
    logger.info(f"Group B bulks: {(bulks_split.obs.target == 'NO-TARGET').sum()}")
    
    # Prepare data
    bulks_split.obs['log10_n_cells'] = np.log10(bulks_split.obs['n_cells'])
    
    # Fit DESeq2 model
    logger.info("Fitting DESeq2 model...")
    model = PyDESeq2(bulks_split, design='~ log10_n_cells + donor + target + ribosomal')
    model.fit(n_cpus=n_cpus_fit)
    
    # Test contrast
    logger.info("Testing contrast: targeting vs NO-TARGET...")
    contrasts = {'targeting_vs_NO-TARGET': (model.cond(target='targeting') - model.cond(target='NO-TARGET'))}
    result_df = model.test_contrasts(contrasts, n_cpus=n_cpus_test)
    
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
    n_sig_005 = (result_df['adj_p_value'] < 0.05).sum()
    n_sig_01 = (result_df['adj_p_value'] < 0.1).sum()
    logger.info(f"Significant genes (p<0.05): {n_sig_005}")
    logger.info(f"Significant genes (p<0.1): {n_sig_01}")
    
    logger.info(f"Iteration {iteration} complete!")

if __name__ == "__main__":
    main()


