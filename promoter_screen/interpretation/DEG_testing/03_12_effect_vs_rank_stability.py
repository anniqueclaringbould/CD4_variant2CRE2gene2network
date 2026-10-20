# %%
"""
One figure for the claim that the estimated effects are stable when a donor is excluded.

Agreement of the log2 fold changes with the full analysis, on the edges reported in the manuscript,
for the three leave-one-donor-out analyses and for the matched random cell subsamples. The random
subsamples are the power-only reference: they remove a comparable number of cells while retaining all
three donors, so a leave-one-donor-out analysis falling below them indicates an effect of donor
composition rather than of reduced cell number.

The axis starts at zero deliberately. All values lie between 0.87 and 1.00 and the point of the figure
is that they are uniformly high, so an expanded axis would exaggerate differences that are not the
message.

Reads the caches from 03_06 (lodo_vs_full_lfc_matrix.parquet) and 03_09
(subsample_vs_full_lfc.parquet), no re-scanning of the chunk files.

Outputs (to --output-dir):
    effect_size_stability.tsv
    effect_size_stability.pdf / .png
"""

import os
import re
import sys
import logging
import argparse

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 8

logging.basicConfig(format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s", stream=sys.stdout)
logger = logging.getLogger('utils')
logger.setLevel(logging.INFO)

FULL = 'full'
FDR = 0.1
DONOR_COLORS = {'A': '#2C7FB8', 'B': '#18974C', 'C': '#E67E22'}
RANDOM_COLORS = {90: '#7A7A7A', 75: '#9E9E9E', 50: '#C4C4C4'}

#: Plotted left to right. Spearman rho and the regression slope are still computed and written to
#: the .tsv -- they say the same thing as the Pearson r and were four bars of repetition on a
#: figure whose point is that every measure agrees.
METRICS = ['Pearson r', 'sign concordance', 'AUCC']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor/output")
    return p.parse_args()


def aucc(px, py):
    """Area under the concordance curve of two adjusted-p rankings of the same edges.

    C(k) is how many of the top k edges the two runs share; AUCC is the mean of C(k)/k over
    k = 1..K. Computed in closed form rather than by looping over k: an edge whose worse rank in
    the two lists is m enters C(k) for every k >= m, so sum_k C(k) = sum_edges (K - m + 1), and
    sum_k k = K(K+1)/2.

    Both lists rank the same K edges here, so this measures whether a leave-one-donor-out analysis
    orders the published edges the way the full analysis does. It is not a recovery rate: an edge
    cannot be missing from a ranking of a fixed edge set.
    """
    ok = np.isfinite(px) & np.isfinite(py)
    K = int(ok.sum())
    mmax = np.maximum(stats.rankdata(px[ok], method='ordinal'),
                      stats.rankdata(py[ok], method='ordinal'))
    return float(np.sum(K - mmax + 1) / (K * (K + 1) / 2))


def agreement(x, y, px, py):
    ok = np.isfinite(x) & np.isfinite(y)
    return {
        'n_edges': int(ok.sum()),
        'Pearson r': float(stats.pearsonr(x[ok], y[ok])[0]),
        'Spearman rho': float(stats.spearmanr(x[ok], y[ok])[0]),
        'sign concordance': float((np.sign(x[ok]) == np.sign(y[ok])).mean()),
        'regression slope': float(np.sum(x[ok] * y[ok]) / np.sum(x[ok] ** 2)),
        'AUCC': aucc(px, py),
    }


def main():
    args = parse_args()
    od = args.output_dir
    f_lodo = os.path.join(od, "lodo_vs_full_lfc_matrix.parquet")
    f_sub = os.path.join(od, "subsample_vs_full_lfc.parquet")
    if not os.path.exists(f_lodo):
        logger.error(f"{f_lodo} not found, run 03_06_compare_lodo_to_full_run.py first.")
        sys.exit(1)

    lfc = pd.read_parquet(f_lodo)
    pub = lfc[lfc[f'sig_{FULL}'].values]
    ref = pub[f'lfc_{FULL}'].values
    ref_p = pub[f'adjp_{FULL}'].values
    logger.info(f"Published edges: {len(pub):,}")

    rows = []
    lodo = sorted(c.replace('lfc_', '') for c in lfc.columns
                  if c.startswith('lfc_') and c != f'lfc_{FULL}')
    for d in lodo:
        rows.append({'run': f'without {d}', 'kind': 'leave one donor out', 'key': d,
                     **agreement(ref, pub[f'lfc_{d}'].values,
                                 ref_p, pub[f'adjp_{d}'].values)})

    if os.path.exists(f_sub):
        sub = pd.read_parquet(f_sub).reindex(pub.index)
        for c in [c for c in sub.columns if c.startswith('lfc_')]:
            pct = int(re.search(r'keep(\d+)', c).group(1))
            rows.append({'run': f'random {pct}%', 'kind': 'random subsample', 'key': pct,
                         **agreement(ref, sub[c].values,
                                     ref_p, sub[c.replace('lfc_', 'adjp_')].values)})
    else:
        logger.warning(f"{f_sub} not found, plotting the leave-one-donor-out analyses only")

    tbl = pd.DataFrame(rows)
    # leave-one-donor-out first, then the random subsamples from least to most cells removed
    tbl['order'] = [0 + i for i in range(len(lodo))] + \
                   [10 - int(k) for k in tbl.loc[tbl.kind == 'random subsample', 'key']]
    tbl = tbl.sort_values(['kind', 'order']).reset_index(drop=True)
    tbl.to_csv(os.path.join(od, "effect_size_stability.tsv"), sep='\t', index=False)
    logger.info("\n" + tbl[['run', 'kind', 'n_edges', 'Pearson r', 'Spearman rho',
                            'sign concordance', 'regression slope', 'AUCC']]
                .round(3).to_string(index=False))

    # %%
    # one wide bar plot: every run x metric as its own bar, laid out left to right in metric blocks
    runs = list(tbl['run'])
    n_r = len(runs)
    gap = 1.0
    pos, labels, colors, vals = [], [], [], []
    block_centres = []
    cur = 0.0
    for m in METRICS:
        start = cur
        for _, r in tbl.iterrows():
            pos.append(cur)
            vals.append(r[m])
            labels.append(r['run'])
            colors.append(DONOR_COLORS.get(str(r['key']), '#999999')
                          if r['kind'] == 'leave one donor out'
                          else RANDOM_COLORS.get(int(r['key']), '#B0B0B0'))
            cur += 1.0
        block_centres.append((start + cur - 1.0) / 2)
        cur += gap

    fig, ax = plt.subplots(figsize=(5.6, 2.4))
    ax.bar(pos, vals, width=0.86, color=colors, edgecolor='white', linewidth=0.4)
    for p, v in zip(pos, vals):
        ax.text(p, v + 0.015, f'{v:.2f}', ha='center', va='bottom', fontsize=5.6, rotation=90)

    ax.set_xticks(pos)
    ax.set_xticklabels(labels, rotation=45, ha='right', rotation_mode='anchor',
                       fontsize=6)
    ax.set_xlim(-0.8, cur - gap - 0.2)
    for c, m in zip(block_centres, METRICS):
        ax.text(c, 1.17, m, ha='center', va='bottom', fontsize=7.5, fontweight='bold')
    # faint separators between the metric blocks
    for c in block_centres[:-1]:
        ax.axvline(c + (n_r - 1) / 2 + gap / 2, color='#DDDDDD', lw=0.9)

    ax.set_ylabel('agreement with the full analysis')
    ax.set_ylim(0, 1.25)
    ax.set_yticks(np.arange(0, 1.01, 0.2))
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    for ext in ['pdf', 'png']:
        out = os.path.join(od, f"effect_size_stability.{ext}")
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Wrote {out}")

    logger.info("All done!")


if __name__ == "__main__":
    main()
