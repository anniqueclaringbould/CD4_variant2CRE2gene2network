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
@click.option('--chunk-start', type=int, required=True, help='Starting chunk index (0-based)')
@click.option('--chunk-end', type=int, default=None, help='Ending chunk index (exclusive). If not provided, processes only chunk-start')
@click.option('--chunk-size', type=int, default=50, help='Number of targets per chunk')
@click.option('--min-cells-per-bulk', type=int, default=5, help='Minimum cells per pseudobulk')
@click.option('--min-bulks', type=int, default=5, help='Minimum pseudobulks per target')
@click.option('--n-cpus-fit', type=int, default=20, help='Number of CPUs for model fitting')
@click.option('--n-cpus-test', type=int, default=20, help='Number of CPUs for contrast testing')
@click.option('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/pseudobulk_deseq2_chunks2", help='Output directory for results')
@click.option('--h5ad-file', type=str, default="/g/stegle/schrod/data/T_Cell/DESeq2_donor_guide_bulks.h5ad", help='Path to input h5ad file')
@click.option('--ribosomal', is_flag=True, help='Correct for ribosomal percentage in design formula')
def main(chunk_start, chunk_end, chunk_size, min_cells_per_bulk, min_bulks, n_cpus_fit, n_cpus_test, output_dir, h5ad_file, ribosomal):
    """
    Run DESeq2 differential expression analysis on specified chunk(s) of targets.
    
    Example usage:
        python 01_02_pseudobulk_DESeq2_chunked_parallel.py --chunk-start 0 --chunk-end 5
        python 01_02_pseudobulk_DESeq2_chunked_parallel.py --chunk-start 10
    """
    
    target_column = "target"
    control_group = "NO-TARGET"
    
    # If chunk_end not specified, process only chunk_start
    if chunk_end is None:
        chunk_end = chunk_start + 1
    
    if ribosomal:
        logger.info("Use ribolomal percentage in design formula!")


    logger.info(f"Processing chunks {chunk_start} to {chunk_end - 1} (chunk size: {chunk_size})")
    logger.info(f"Parameters: min_cells_per_bulk={min_cells_per_bulk}, min_bulks={min_bulks}")
    logger.info(f"CPUs: fit={n_cpus_fit}, test={n_cpus_test}")
    
    # %%
    logger.info(f"Read Bulked anndata object from: {h5ad_file}")
    bulks = ad.read_h5ad(os.path.join(h5ad_file))
    logger.info(f"Anndata object with {bulks.n_obs} cells and {bulks.n_vars} genes loaded.")

    # %%
    #remove small bulks
    logger.info(f"Filtering pseudobulks with fewer than {min_cells_per_bulk} cells...")
    n_before = bulks.n_obs
    bulks = bulks[bulks.obs.n_cells>=min_cells_per_bulk]
    logger.info(f"Removed {n_before - bulks.n_obs} pseudobulks. {bulks.n_obs} remaining.")

    # %%
    logger.info(f"Extracting target gene names from guide IDs...")
    bulks.obs["target"] = bulks.obs.guide.apply(lambda x: "-".join(x.split("-")[1:-1]))
    logger.info(f"Found {bulks.obs.target.nunique()} unique targets.")

    # %%
    # Only keep targets with at least min_bulks bulks
    logger.info(f"Filtering targets with at least {min_bulks} pseudobulks...")
    n_before = bulks.n_obs
    valid_targets = bulks.obs.target.value_counts()>=min_bulks
    valid_targets = valid_targets[valid_targets].index.tolist()
    bulks = bulks[bulks.obs.target.isin(valid_targets + ["NO-TARGET"])]
    logger.info(f"Kept {len(valid_targets)} targets with sufficient pseudobulks. {bulks.n_obs} pseudobulks remaining.")

    #%%
    # Run deseq in chunks of targets to avoid memory issues
    all_targets = [t for t in bulks.obs.target.unique() if t != "NO-TARGET"]
    n_chunks = int(np.ceil(len(all_targets) / chunk_size))
    
    logger.info(f"Total targets: {len(all_targets)}, Total chunks: {n_chunks}")
    
    # Validate chunk range
    if chunk_start < 0 or chunk_start >= n_chunks:
        logger.error(f"Invalid chunk_start: {chunk_start}. Must be in range [0, {n_chunks - 1}]")
        sys.exit(1)
    
    if chunk_end > n_chunks:
        logger.warning(f"chunk_end {chunk_end} exceeds total chunks {n_chunks}. Setting to {n_chunks}")
        chunk_end = n_chunks
    
    os.makedirs(output_dir, exist_ok=True)
    logger.info(f"Saving chunk results to {output_dir}/")

    for chunk_idx in range(chunk_start, chunk_end):
        start_idx = chunk_idx * chunk_size
        end_idx = min((chunk_idx + 1) * chunk_size, len(all_targets))
        chunk_targets = all_targets[start_idx:end_idx]
        
        logger.info(f"Chunk {chunk_idx + 1}/{n_chunks}: Processing targets {start_idx + 1}-{end_idx} ({len(chunk_targets)} targets)...")
        
        # Select bulks for this chunk (always include NO-TARGET)
        chunk_bulks = bulks[bulks.obs.target.isin(chunk_targets + ["NO-TARGET"])].copy()
        logger.info(f"  Selected {chunk_bulks.n_obs} pseudobulks ({chunk_bulks.obs.target.value_counts()['NO-TARGET']} NO-TARGET)")
        
        # Prepare data
        chunk_bulks.obs['log10_n_cells'] = np.log10(chunk_bulks.obs['n_cells'])
        
        # Fit model
        logger.info(f"  Fitting DESeq2 model with {chunk_bulks.n_obs} pseudobulks and {chunk_bulks.n_vars} genes...")
        
        if ribosomal:
            model = PyDESeq2(chunk_bulks, design='~ log10_n_cells + donor + target + ribosomal')
        else:
            model = PyDESeq2(chunk_bulks, design='~ log10_n_cells + donor + target')
        
        model.fit(n_cpus=n_cpus_fit)
        
        # Test contrasts
        logger.info(f"  Testing {len(chunk_targets)} contrasts...")
        contrasts = {t:(model.cond(target=t) - model.cond(target="NO-TARGET")) for t in chunk_targets}
        chunk_res_df = model.test_contrasts(contrasts, n_cpus=n_cpus_test)
        logger.info(f"  Chunk results shape: {chunk_res_df.shape}")
        
        # Save chunk results
        chunk_output_path = os.path.join(output_dir, f"chunk_{chunk_idx + 1:03d}_of_{n_chunks:03d}.csv")
        chunk_res_df.to_csv(chunk_output_path)
        logger.info(f"  Saved chunk {chunk_idx + 1} results to {chunk_output_path}")

    logger.info(f"Completed processing chunks {chunk_start} to {chunk_end - 1}")

if __name__ == "__main__":
    main()
