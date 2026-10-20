# %%
"""
Summarise the random cell subsampling runs of 03_01 (--keep-frac), so they can be plotted next to the
leave-one-donor-out (LODO) runs.

The point of the comparison: the LODO runs retain 79.4%, 65.0% and 55.6% of the cells, so part of what
they lose is simply power. Randomly removing a comparable fraction of cells, with all three donors kept,
gives the power-only reference. If the LODO runs fall below the random subsampling curve, the shortfall
is attributable to donor identity rather than to cell number.

Edges are counted on exactly the same universe as the LODO numbers in 03_05: perturbations and
genes tested in all three LODO runs, self edges excluded, adj_p < 0.1. Cell counts are the nominal
keep_frac times the total, which the binomial draw matches to within a few tenths of a percent.

Outputs (to --output-dir):
    subsample_edges_near_significant.parquet   adj_p < 0.5 per subsample run (intermediate)
    subsample_summary.tsv                      keep_frac, cells retained, perturbations, edges
"""

import os
import re
import sys
import argparse

import numpy as np
import pandas as pd

os.chdir("/g/stegle/schrod/code/TCell")
from utils import *

logger = get_logger(__name__)

FDR = 0.1
CACHE_P = 0.5


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--subsample-dir', type=str, default="/g/stegle/schrod/code/TCell/Subsample",
                   help='Directory holding the pseudobulk_deseq2_keep*_seed* result directories')
    p.add_argument('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor/output",
                   help='Where the LODO tables live and where the summary is written')
    p.add_argument('--force-reload', action='store_true')
    return p.parse_args()


def main():
    args = parse_args()
    od = args.output_dir
    cache = os.path.join(od, "subsample_edges_near_significant.parquet")
    tc_f = os.path.join(od, "lodo_tested_contrasts.csv")
    tg_f = os.path.join(od, "lodo_tested_genes.csv")
    dc_f = os.path.join(od, "donor_contribution.tsv")
    for f in [tc_f, tg_f, dc_f]:
        if not os.path.exists(f):
            logger.error(f"{f} not found, run 03_05 and 03_07 first.")
            sys.exit(1)

    # the universe used for every LODO number
    tc, tg = pd.read_csv(tc_f), pd.read_csv(tg_f)
    lodo_cols = [c for c in tc.columns if c != 'contrast']
    common_c = set(tc.loc[tc[lodo_cols].all(axis=1), 'contrast'])
    common_g = set(tg.loc[tg[lodo_cols].all(axis=1), 'gene'])
    total_cells = int(pd.read_csv(dc_f, sep='\t')['total_cells'].sum())
    logger.info(f"Universe: {len(common_c):,} perturbations x {len(common_g):,} genes; "
                f"{total_cells:,} cells in the full data")

    # %%
    runs = sorted(d for d in os.listdir(args.subsample_dir)
                  if d.startswith('pseudobulk_deseq2_keep'))
    if not runs:
        logger.error(f"No pseudobulk_deseq2_keep* directories found in {args.subsample_dir}")
        sys.exit(1)
    logger.info(f"Found subsample runs: {', '.join(runs)}")

    if os.path.exists(cache) and not args.force_reload:
        logger.info(f"Loading cached edges from {cache}")
        near = pd.read_parquet(cache)
    else:
        parts = []
        for d in runs:
            m = re.match(r'pseudobulk_deseq2_keep(\d+)_seed(\d+)', d)
            keep_pct, seed = int(m.group(1)), int(m.group(2))
            path = os.path.join(args.subsample_dir, d)
            files = sorted(f for f in os.listdir(path) if f.endswith('.csv'))
            logger.info(f"  {d}: reading {len(files)} chunk files")
            keep = []
            for i, fn in enumerate(files, start=1):
                df = pd.read_csv(os.path.join(path, fn),
                                 usecols=['contrast', 'variable', 'adj_p_value'])
                keep.append(df[df['adj_p_value'] < CACHE_P])
                if i % 50 == 0 or i == len(files):
                    logger.info(f"    {i}/{len(files)} files")
            e = pd.concat(keep, ignore_index=True)
            e['keep_pct'], e['seed'] = keep_pct, seed
            parts.append(e)
        near = pd.concat(parts, ignore_index=True)
        near['edge'] = near['contrast'] + '->' + near['variable']
        near = near[['edge', 'adj_p_value', 'keep_pct', 'seed']]
        near.to_parquet(cache, index=False)
        logger.info(f"Cached {len(near):,} rows (adj_p < {CACHE_P}) to {cache}")

    # %%
    # restrict to the LODO universe and count edges
    near = near.copy()
    near['contrast'] = near['edge'].str.split('->').str[0]
    near['variable'] = near['edge'].str.split('->').str[1]
    near = near[near['contrast'].isin(common_c) & near['variable'].isin(common_g)]
    near = near[near['contrast'] != near['variable']]

    rows = []
    for (kp, seed), g in near.groupby(['keep_pct', 'seed']):
        sig = g[g['adj_p_value'] < FDR]
        keep_frac = kp / 100
        rows.append({'keep_pct': kp, 'seed': seed, 'keep_frac': keep_frac,
                     'cells_retained': int(round(keep_frac * total_cells)),
                     'pct_cells_retained': round(100 * keep_frac, 1),
                     'perturbations_with_edge': sig['contrast'].nunique(),
                     'edges_detected': len(set(sig['edge']))})
    summary = pd.DataFrame(rows).sort_values('keep_pct')
    summary.to_csv(os.path.join(od, "subsample_summary.tsv"), sep='\t', index=False)
    logger.info("\n" + summary.to_string(index=False))
    logger.info("All done!")


if __name__ == "__main__":
    main()
