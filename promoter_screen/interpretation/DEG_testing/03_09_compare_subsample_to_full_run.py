# %%
"""
Same comparison as 03_06, but for the random cell subsampling runs, plotted next to the
leave-one-donor-out (LODO) runs so the two kinds of data loss can be read off the same axes.

The LODO runs retain 79.4%, 65.0% and 55.6% of the cells. The random subsamples retain 90%, 75% and 50%
with all three donors kept, so they are the power-only reference: any difference between a LODO run and
the random curve at a matched cell count is attributable to which cells were removed rather than how
many.

Recall is the fraction of published edges a run recovers, precision the fraction of a run's calls that
are published edges. Both are computed on the universe used throughout: perturbations and genes tested
in all three LODO runs, self edges excluded, FDR 0.1.

Requires the caches from 03_06 (lodo_vs_full_lfc_matrix.parquet), 03_05
(all_runs_near_significant.parquet) and 03_08 (subsample_edges_near_significant.parquet).

Outputs (to --output-dir):
    subsample_vs_full_lfc.parquet          LFC of the published edges in each subsample run
    subsample_vs_full_recovery.tsv         recall, precision and LFC agreement for every run
    subsample_vs_full_comparison.pdf/.png
"""

import os
import re
import sys
import argparse
import itertools

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 8

os.chdir("/g/stegle/schrod/code/TCell")
from utils import *

logger = get_logger(__name__)

FULL = 'full'
FDR = 0.1
C_RECALL, C_PREC = '#2C7FB8', '#18974C'
C_LODO, C_SUB = '#D62728', '#555555'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--subsample-dir', type=str, default="/g/stegle/schrod/code/TCell/Subsample")
    p.add_argument('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor/output")
    p.add_argument('--force-reload', action='store_true')
    return p.parse_args()


def main():
    args = parse_args()
    od = args.output_dir
    lfc_f = os.path.join(od, "lodo_vs_full_lfc_matrix.parquet")
    near_f = os.path.join(od, "all_runs_near_significant.parquet")
    sub_near_f = os.path.join(od, "subsample_edges_near_significant.parquet")
    tc_f, tg_f = (os.path.join(od, "lodo_tested_contrasts.csv"),
                  os.path.join(od, "lodo_tested_genes.csv"))
    for f in [lfc_f, near_f, sub_near_f, tc_f, tg_f]:
        if not os.path.exists(f):
            logger.error(f"{f} not found, run 03_05 / 03_06 / 03_08 first.")
            sys.exit(1)
    cache = os.path.join(od, "subsample_vs_full_lfc.parquet")

    # universe
    tc, tg = pd.read_csv(tc_f), pd.read_csv(tg_f)
    lodo = [c for c in tc.columns if c != 'contrast']
    common_c = set(tc.loc[tc[lodo].all(axis=1), 'contrast'])
    common_g = set(tg.loc[tg[lodo].all(axis=1), 'gene'])

    # published edge set and the LODO LFCs
    lfc = pd.read_parquet(lfc_f)
    pub = lfc[lfc[f'sig_{FULL}'].values]
    published = set(pub.index)
    logger.info(f"Published edges: {len(published):,}")

    # %%
    # Pass over the subsample runs for the published edges' LFC and adj_p
    runs_sub = sorted(d for d in os.listdir(args.subsample_dir)
                      if d.startswith('pseudobulk_deseq2_keep'))
    keep_of = {d: int(re.match(r'pseudobulk_deseq2_keep(\d+)_seed(\d+)', d).group(1)) for d in runs_sub}
    sub_labels = {d: f'random keep {keep_of[d]}%' for d in runs_sub}

    if os.path.exists(cache) and not args.force_reload:
        logger.info(f"Loading cached subsample LFCs from {cache}")
        sub_wide = pd.read_parquet(cache)
    else:
        cols = []
        for d in runs_sub:
            path = os.path.join(args.subsample_dir, d)
            files = sorted(f for f in os.listdir(path) if f.endswith('.csv'))
            logger.info(f"  {d}: extracting {len(published):,} published edges from {len(files)} files")
            out = []
            for i, fn in enumerate(files, start=1):
                df = pd.read_csv(os.path.join(path, fn),
                                 usecols=['contrast', 'variable', 'log_fc', 'adj_p_value'])
                df['edge'] = df['contrast'] + '->' + df['variable']
                out.append(df.loc[df['edge'].isin(published), ['edge', 'log_fc', 'adj_p_value']])
                if i % 50 == 0 or i == len(files):
                    logger.info(f"    {i}/{len(files)} files")
            r = pd.concat(out, ignore_index=True).set_index('edge')
            cols.append(r.rename(columns={'log_fc': f'lfc_{d}', 'adj_p_value': f'adjp_{d}'}))
        sub_wide = pd.concat(cols, axis=1)
        sub_wide.to_parquet(cache)
        logger.info(f"Cached {sub_wide.shape[0]:,} published edges x {len(runs_sub)} runs to {cache}")

    # %%
    # Edge sets at FDR 0.1 for every run, on the common universe
    def restrict(df, run_col, run):
        d = df[df[run_col] == run].copy()
        d['contrast'] = d['edge'].str.split('->').str[0]
        d['variable'] = d['edge'].str.split('->').str[1]
        d = d[d['contrast'].isin(common_c) & d['variable'].isin(common_g)]
        d = d[d['contrast'] != d['variable']]
        return set(d.loc[d['adj_p_value'] < FDR, 'edge'])

    near = pd.read_parquet(near_f)
    E = {r: restrict(near, 'run', r) for r in [FULL] + lodo}
    sub_near = pd.read_parquet(sub_near_f)
    for d in runs_sub:
        E[d] = restrict(sub_near, 'keep_pct', keep_of[d])

    # %%
    # Recall, precision and LFC agreement per run
    order = lodo + runs_sub
    labels = {**{r: f'without {r}' for r in lodo}, **sub_labels}
    rows = []
    for r in order:
        shared = len(E[r] & E[FULL])
        col = f'lfc_{r}' if r in lodo else f'lfc_{r}'
        x = pub[f'lfc_{FULL}'].values
        y = (pub[col].values if r in lodo else
             sub_wide.reindex(pub.index)[f'lfc_{r}'].values)
        ok = np.isfinite(x) & np.isfinite(y)
        rows.append({
            'run': labels[r], 'kind': 'leave one donor out' if r in lodo else 'random subsample',
            'n_edges': len(E[r]), 'shared_with_published': shared,
            'recall_pct': round(100 * shared / len(E[FULL]), 1),
            'precision_pct': round(100 * shared / len(E[r]), 1) if E[r] else np.nan,
            'jaccard': round(shared / len(E[r] | E[FULL]), 3),
            'lfc_pearson_r': round(float(stats.pearsonr(x[ok], y[ok])[0]), 3),
            'lfc_spearman_rho': round(float(stats.spearmanr(x[ok], y[ok])[0]), 3),
            'lfc_sign_concordance_pct': round(100 * float((np.sign(x[ok]) == np.sign(y[ok])).mean()), 1),
        })
    tbl = pd.DataFrame(rows)
    tbl.to_csv(os.path.join(od, "subsample_vs_full_recovery.tsv"), sep='\t', index=False)
    logger.info("\n" + tbl.to_string(index=False))

    # %%
    fig = plt.figure(figsize=(11.5, 7.6))
    outer = gridspec.GridSpec(2, 2, figure=fig, hspace=0.55, wspace=0.3, height_ratios=[1, 1])

    # A, recall and precision for every run
    ax = fig.add_subplot(outer[0, 0])
    x = np.arange(len(tbl))
    ax.bar(x - 0.2, tbl['recall_pct'], width=0.4, color=C_RECALL,
           label='recall: % of published edges recovered')
    ax.bar(x + 0.2, tbl['precision_pct'], width=0.4, color=C_PREC,
           label='precision: % of calls that are published edges')
    for xi, (a, b) in enumerate(zip(tbl['recall_pct'], tbl['precision_pct'])):
        ax.text(xi - 0.2, a, f'{a:.0f}', ha='center', va='bottom', fontsize=6.5)
        ax.text(xi + 0.2, b, f'{b:.0f}', ha='center', va='bottom', fontsize=6.5)
    ax.set_xticks(x)
    ax.set_xticklabels(tbl['run'], rotation=30, ha='right', fontsize=7)
    ax.axvline(len(lodo) - 0.5, color='#999999', ls=':', lw=1)
    ax.set_ylabel('Percent')
    ax.set_ylim(0, 108)
    ax.legend(fontsize=6.5, frameon=False, loc='upper center')
    ax.set_title(f'A  Recall and precision vs the published run (FDR < {FDR})',
                 fontsize=9, fontweight='bold', loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    # B, LFC Pearson r for every run
    ax = fig.add_subplot(outer[0, 1])
    ax.bar(x, tbl['lfc_pearson_r'], width=0.6,
           color=[C_LODO if k == 'leave one donor out' else C_SUB for k in tbl['kind']])
    for xi, v, s in zip(x, tbl['lfc_pearson_r'], tbl['lfc_sign_concordance_pct']):
        ax.text(xi, v, f'{v:.2f}\nsign {s:.0f}%', ha='center', va='bottom', fontsize=6.5)
    ax.set_xticks(x)
    ax.set_xticklabels(tbl['run'], rotation=30, ha='right', fontsize=7)
    ax.axvline(len(lodo) - 0.5, color='#999999', ls=':', lw=1)
    ax.set_ylabel('LFC Pearson r vs the full run')
    ax.set_ylim(0, 1.18)
    ax.set_title('B  Effect size agreement on the published edges',
                 fontsize=9, fontweight='bold', loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    # C, LFC hexbins for the random subsamples
    inner = gridspec.GridSpecFromSubplotSpec(1, len(runs_sub), subplot_spec=outer[1, :], wspace=0.3)
    for i, d in enumerate(runs_sub):
        ax = fig.add_subplot(inner[i])
        xv = pub[f'lfc_{FULL}'].values
        yv = sub_wide.reindex(pub.index)[f'lfc_{d}'].values
        ok = np.isfinite(xv) & np.isfinite(yv)
        ax.hexbin(xv[ok], yv[ok], gridsize=55, bins='log', cmap='Blues', mincnt=1)
        lim = np.nanpercentile(np.abs(np.concatenate([xv[ok], yv[ok]])), 99.5)
        ax.plot([-lim, lim], [-lim, lim], color='#D62728', ls='--', lw=1)
        ax.axhline(0, color='#B0B0B0', lw=0.5)
        ax.axvline(0, color='#B0B0B0', lw=0.5)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        rr = stats.pearsonr(xv[ok], yv[ok])[0]
        rho = stats.spearmanr(xv[ok], yv[ok])[0]
        conc = 100 * (np.sign(xv[ok]) == np.sign(yv[ok])).mean()
        ax.set_xlabel('LFC full run (3 donors)')
        ax.set_ylabel(f'LFC {sub_labels[d]}' if i == 0 else '')
        ax.set_title(f'{"C  " if i == 0 else ""}{sub_labels[d]}\n'
                     f'r={rr:.2f}, rho={rho:.2f}, sign {conc:.0f}%  (n={ok.sum():,})',
                     fontsize=8, fontweight='bold' if i == 0 else 'normal', loc='left')
        for s in ['top', 'right']:
            ax.spines[s].set_visible(False)

    for ext in ['pdf', 'png']:
        out = os.path.join(od, f"subsample_vs_full_comparison.{ext}")
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Wrote {out}")

    logger.info("All done!")


if __name__ == "__main__":
    main()
