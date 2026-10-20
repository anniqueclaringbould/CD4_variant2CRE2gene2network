# %%
"""
Aggregate the per-cell AUC scores to the pseudobulk level, so they can be used as covariates in
03_01_run_DESeq2_batched.py.

The scores are the cell state signatures already used in the single cell limma design
(CellCycle.G2M, early, ICOS.CD38, Mito, CTLA4.CD38, Poor.Quality, iTreg_5d4h, Th1_5d). Each pseudobulk
gets the mean score of the cells that went into it, which matches how the existing `ribosomal` covariate
is defined.

Written as a side table rather than into the pseudobulk h5ad, so that the published pseudobulks are not
modified and the leave-one-donor-out and subsampled runs can use the same file without rebuilding
anything. Only .obs is read from the 56 GB cell level h5ad, never the counts.

Score names are sanitised (dots to underscores: CellCycle.G2M -> CellCycle_G2M) because a dot in a bare
variable name breaks formula parsing in the DESeq2 design.

Outputs:
    --out-file (default /g/stegle/schrod/data/T_Cell/pseudobulk_cell_scores.tsv)
        one row per pseudobulk, indexed to match the pseudobulk .obs index (e.g. A_GUIDE-A1BG-2)
"""

import os
import sys
import argparse

import h5py
import numpy as np
import pandas as pd

os.chdir("/g/stegle/schrod/code/TCell")
from utils import *

logger = get_logger(__name__)

SCORE_FILE = ("/g/steinmetz/project/otar_2063/promoter_screen/interpretation/cell_scoring/"
              "output/Representative_AUC_scores_per_cell.tsv")
CELL_H5AD = "/g/stegle/schrod/data/T_Cell/promoter_single_guide_hvg_raw_counts.h5ad"
BULK_H5AD = "/g/stegle/schrod/data/T_Cell/DESeq2_donor_guide_bulks.h5ad"
OUT_FILE = "/g/stegle/schrod/data/T_Cell/pseudobulk_cell_scores.tsv"


def read_obs(path, cols):
    """Read selected .obs columns (and the index) from an h5ad without touching X."""
    out = {}
    with h5py.File(path, 'r') as f:
        idx_key = f['obs'].attrs.get('_index', '_index')
        index = pd.Index([x.decode() for x in f['obs'][idx_key][:]])
        for c in cols:
            g = f['obs'][c]
            if isinstance(g, h5py.Group):        # categorical
                cats = g['categories'][:].astype(str)
                out[c] = pd.Categorical.from_codes(g['codes'][:], cats).astype(str)
            else:
                out[c] = g[:]
    return pd.DataFrame(out, index=index)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--score-file', type=str, default=SCORE_FILE)
    p.add_argument('--cell-h5ad', type=str, default=CELL_H5AD)
    p.add_argument('--bulk-h5ad', type=str, default=BULK_H5AD,
                   help='Only used to check that the produced index matches the pseudobulks')
    p.add_argument('--out-file', type=str, default=OUT_FILE)
    return p.parse_args()


def main():
    args = parse_args()

    # %%
    logger.info(f"Reading per-cell scores from {args.score_file}")
    scores = pd.read_csv(args.score_file, sep='\t').set_index('cell')
    logger.info(f"  {scores.shape[0]:,} cells x {scores.shape[1]} scores: "
                f"{', '.join(scores.columns)}")

    # %%
    logger.info(f"Reading cell metadata from {args.cell_h5ad}")
    obs = read_obs(args.cell_h5ad, ['donor', 'guide'])
    logger.info(f"  {len(obs):,} cells")

    missing = (~obs.index.isin(scores.index)).sum()
    if missing:
        logger.warning(f"{missing:,} of {len(obs):,} cells have no score and are dropped from the "
                       f"pseudobulk means")
    obs = obs.join(scores, how='inner')
    logger.info(f"  {len(obs):,} cells with a score")

    # %%
    # the same key 01_02_build_pseudobulks.py uses to define a pseudobulk
    obs['sample_id'] = obs['donor'].astype(str) + '_' + obs['guide'].astype(str)

    logger.info("Averaging scores per pseudobulk...")
    agg = obs.groupby('sample_id')[list(scores.columns)].mean()
    agg.insert(0, 'n_cells_scored', obs.groupby('sample_id').size())

    # dots break formula parsing in the design, so make the names formula safe
    renamed = {c: c.replace('.', '_') for c in scores.columns}
    agg = agg.rename(columns=renamed)
    changed = {k: v for k, v in renamed.items() if k != v}
    if changed:
        logger.info(f"  renamed for formula safety: {changed}")

    # %%
    # sanity check against the pseudobulks the scores will be joined to
    if os.path.exists(args.bulk_h5ad):
        with h5py.File(args.bulk_h5ad, 'r') as f:
            idx_key = f['obs'].attrs.get('_index', '_index')
            bulk_idx = pd.Index([x.decode() for x in f['obs'][idx_key][:]])
        n_hit = bulk_idx.isin(agg.index).sum()
        logger.info(f"Pseudobulks in {os.path.basename(args.bulk_h5ad)}: {len(bulk_idx):,}, "
                    f"of which {n_hit:,} ({100 * n_hit / len(bulk_idx):.2f}%) have scores")
        if n_hit < len(bulk_idx):
            logger.warning("Some pseudobulks have no scores, they would become NaN in the design")

    agg.index.name = 'sample_id'
    agg.to_csv(args.out_file, sep='\t')
    logger.info(f"Wrote {agg.shape[0]:,} pseudobulks x {agg.shape[1] - 1} scores to {args.out_file}")
    logger.info("\n" + agg.describe().T.round(3).to_string())
    logger.info("All done!")


if __name__ == "__main__":
    main()
