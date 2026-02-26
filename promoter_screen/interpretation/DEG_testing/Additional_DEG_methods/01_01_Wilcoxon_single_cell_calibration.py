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
import itertools
from tqdm import tqdm

from utils import *

import logging

logging.basicConfig(format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s", stream=sys.stdout)
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)



#%%
# Run if NT-adata does not exist
subset = False
if subset:
    adata_path = "/g/stegle/schrod/data/T_Cell/promoter_single_guide_hvg.h5ad"
    adata = ad.read_h5ad(os.path.join(adata_path))
    adata_NT = adata[adata.obs.target=="NO-TARGET"]
    adata_NT.write_h5ad("/g/stegle/schrod/data/T_Cell/promoter_NT_hvg.h5ad")

#%%
min_cells = 5
target_column = "guide"
h5ad_file = f"/g/stegle/schrod/data/T_Cell/promoter_NT_hvg.h5ad"

#%%
logger.info(f"Read anndata object from: {h5ad_file}")
adata_NT = ad.read_h5ad(os.path.join(h5ad_file))
adata_NT.obs["donor"] = adata_NT.obs.batch.apply(lambda x: x[0])

#%%
counts = adata_NT.obs["guide"].value_counts()
selected_guides = counts[counts > min_cells].index
adata_NT = adata_NT[adata_NT.obs.guide.isin(selected_guides)].copy()
guides = adata_NT.obs.guide.value_counts()


#%%
logger.info(f"Computing differentially expressed genes...")
sc.tl.rank_genes_groups(adata_NT, groupby=target_column, groups=list(guides.index), reference="rest", method='wilcoxon')
logger.info(f"Finished computing differentially expressed genes!")
df = get_diff_expression_results(adata_NT)
df.to_csv(f"/g/stegle/schrod/code/TCell/results/NT_one_against_rest_DEG_results.csv")