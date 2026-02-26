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

os.chdir("/g/stegle/schrod/code/TCell")
from utils import *

logger = get_logger(__name__)

# %%
min_cells_per_bulk = 5
target_column = "target"
control_group = "NO-TARGET"
h5ad_file = f"/g/stegle/schrod/data/T_Cell/promoter_single_guide_hvg_raw_counts.h5ad"

# %%
logger.info(f"Read anndata object from: {h5ad_file}")
adata = ad.read_h5ad(os.path.join(h5ad_file))
logger.info(f"Anndata object with {adata.n_obs} cells and {adata.n_vars} genes loaded.")

# %%
logger.info(f"Build pseudobulk profiles by aggregating cells with the same guide per donor...")
sample_cols = ['donor', 'guide']   # Columns used to aggregate cells
adata.obs["sample_id"] = adata.obs[sample_cols].apply(lambda x: "_".join(x), axis=1)



# %%
logger.info(f"Aggregating cells into pseudobulk profiles...")
n_cells_obs = adata.obs.value_counts(['sample_id'] + sample_cols).reset_index()
n_cells_obs = n_cells_obs.set_index('sample_id').rename({'count':'n_cells'}, axis=1)
bulks = sc.get.aggregate(adata, by='sample_id', func='sum')
bulks.obs = n_cells_obs.loc[bulks.obs_names].copy()
logger.info(f"Created {bulks.n_obs} pseudobulk profiles.")

#%%
bulks.X = bulks.layers['sum']

# %%
# Get genes starting with RPL and RPS
gene_names = bulks.var_names.to_list()
rpl_genes = [gene for gene in gene_names if gene.startswith("RPL")]
rps_genes = [gene for gene in gene_names if gene.startswith("RPS")]
print(f"Found {len(rpl_genes)} RPL genes and {len(rps_genes)} RPS genes.")

# For each bulk get all sum all of these genes
bulks.obs["ribosomal"] = (bulks[:, rpl_genes].X.sum(axis=1) + bulks[:, rps_genes].X.sum(axis=1))/bulks.X.sum(1)

#%%
# Save the bulks anndata object for downstream analysis
bulks.write_h5ad("/g/stegle/schrod/data/T_Cell/DESeq2_donor_guide_bulks.h5ad")
# %%
