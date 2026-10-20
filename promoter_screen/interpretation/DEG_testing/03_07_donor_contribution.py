# %%
"""
Why leaving out one donor changes the results as much as it does.

The leave-one-donor-out (LODO) runs differ far more than the "three donors" framing suggests: leaving
out donor A reproduces 98% of the full run's edge count while leaving out donor C reproduces only 48%.
This script shows that the donors contribute very unequal amounts of data, and relates each LODO run's
discovery count to the fraction of cells it retains.

Note that the *number* of pseudobulks is nearly identical per donor by construction (3 gRNAs x ~7,000
targets), so it is the cells per pseudobulk, not the pseudobulk count, that differs.

Reads the pseudobulk AnnData that DESeq2 was run on, via h5py so the 7.4 GB counts matrix is never
loaded.

Outputs (to --output-dir):
    donor_contribution.tsv          per donor cell and pseudobulk statistics
    donor_contribution.pdf / .png
"""

import os
import argparse

import h5py
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 8

os.chdir("/g/stegle/schrod/code/TCell")
from utils import *

logger = get_logger(__name__)

MIN_CELLS = 5          # the pseudobulk filter used in 03_01
# edges detected at FDR 0.1 on the commonly tested perturbations, from 03_05
EDGES = {'full': 74296, 'A': 72813, 'B': 47282, 'C': 35804}
TESTED = {'full': 7300, 'A': 6778, 'B': 6623, 'C': 6146}

COLORS = {'A': '#2C7FB8', 'B': '#18974C', 'C': '#E67E22'}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--h5ad', type=str,
                   default="/g/stegle/schrod/data/T_Cell/DESeq2_donor_guide_bulks.h5ad")
    p.add_argument('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor/output")
    return p.parse_args()


def read_obs(path, cols):
    """Read selected .obs columns from an h5ad without touching X."""
    out = {}
    with h5py.File(path, 'r') as f:
        for c in cols:
            g = f['obs'][c]
            if isinstance(g, h5py.Group):      # categorical
                cats = g['categories'][:].astype(str)
                out[c] = pd.Categorical.from_codes(g['codes'][:], cats)
            else:
                out[c] = g[:]
    return pd.DataFrame(out)


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    obs = read_obs(args.h5ad, ['donor', 'guide', 'n_cells', 'ribosomal'])
    logger.info(f"{len(obs):,} pseudobulks (donor x gRNA)")

    g = obs.groupby('donor', observed=True)
    tbl = pd.DataFrame({
        'n_pseudobulks': g.size(),
        'total_cells': g['n_cells'].sum(),
        'median_cells_per_pseudobulk': g['n_cells'].median(),
        'mean_cells_per_pseudobulk': g['n_cells'].mean().round(1),
        f'pseudobulks_below_{MIN_CELLS}_cells': g['n_cells'].apply(lambda v: int((v < MIN_CELLS).sum())),
        f'pct_dropped_by_min_cells': (g['n_cells'].apply(lambda v: (v < MIN_CELLS).mean()) * 100).round(1),
        'median_ribosomal_fraction': g['ribosomal'].median().round(4),
    })
    tbl['pct_of_all_cells'] = (100 * tbl['total_cells'] / tbl['total_cells'].sum()).round(1)

    # what each LODO run retains, and what it found
    total_cells = tbl['total_cells'].sum()
    rows = []
    for d in tbl.index:
        kept = total_cells - tbl.loc[d, 'total_cells']
        rows.append({'left_out_donor': d,
                     'cells_retained': int(kept),
                     'pct_cells_retained': round(100 * kept / total_cells, 1),
                     'perturbations_tested': TESTED[d],
                     'edges_detected': EDGES[d],
                     'pct_of_full_run_edges': round(100 * EDGES[d] / EDGES['full'], 1)})
    lodo = pd.DataFrame(rows)

    tbl.to_csv(os.path.join(args.output_dir, "donor_contribution.tsv"), sep='\t')
    logger.info("\nPer donor:\n" + tbl.to_string())
    logger.info("\nPer LODO run:\n" + lodo.to_string(index=False))

    # random cell subsampling reference, if 03_08_subsample_summary.py has been run. All three donors are kept, so this is
    # the power-only curve that the LODO runs should be compared against.
    sub_f = os.path.join(args.output_dir, "subsample_summary.tsv")
    sub = pd.read_csv(sub_f, sep='\t') if os.path.exists(sub_f) else None
    if sub is not None:
        sub = pd.concat([sub, pd.DataFrame([{'pct_cells_retained': 100.0,
                                             'edges_detected': EDGES['full']}])],
                        ignore_index=True).sort_values('pct_cells_retained')
        logger.info("\nRandom subsampling reference:\n" +
                    sub[['pct_cells_retained', 'edges_detected']].to_string(index=False))
        # how far each LODO run sits from the random curve at the same cell count
        for _, r in lodo.iterrows():
            exp = np.interp(r['pct_cells_retained'], sub['pct_cells_retained'], sub['edges_detected'])
            logger.info(f"  without {r['left_out_donor']}: {r['edges_detected']:,} edges vs "
                        f"{exp:,.0f} expected from random subsampling at "
                        f"{r['pct_cells_retained']:.1f}% cells "
                        f"({100 * (r['edges_detected'] / exp - 1):+.1f}%)")
    else:
        logger.warning(f"{sub_f} not found, plotting without the random subsampling reference")

    # %%
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4))

    # A, cells per donor
    ax = axes[0]
    donors = list(tbl.index)
    ax.bar(donors, tbl['total_cells'] / 1e3, color=[COLORS[d] for d in donors], width=0.6)
    for i, d in enumerate(donors):
        ax.text(i, tbl.loc[d, 'total_cells'] / 1e3,
                f"{tbl.loc[d, 'pct_of_all_cells']:.0f}%", ha='center', va='bottom', fontsize=8)
    ax.set_ylabel('Cells in pseudobulks (thousands)')
    ax.set_xlabel('Donor')
    ax.set_title('A  Cells contributed per donor', fontsize=9, fontweight='bold', loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    # B, cells per pseudobulk distribution
    ax = axes[1]
    data = [obs.loc[obs['donor'] == d, 'n_cells'].values for d in donors]
    bp = ax.boxplot(data, labels=donors, showfliers=False, widths=0.55, patch_artist=True,
                    medianprops=dict(color='#111111'))
    for patch, d in zip(bp['boxes'], donors):
        patch.set_facecolor(COLORS[d])
        patch.set_alpha(0.65)
    ax.axhline(MIN_CELLS, color='#D62728', ls='--', lw=1.1)
    ax.text(0.55, MIN_CELLS, f' min {MIN_CELLS} cells', color='#D62728', fontsize=7, va='bottom')
    ax.set_yscale('log')
    ax.set_ylabel('Cells per pseudobulk')
    ax.set_xlabel('Donor')
    ax.set_title('B  Pseudobulk depth per donor', fontsize=9, fontweight='bold', loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    # C, discoveries vs cells retained, LODO against the random subsampling reference
    ax = axes[2]
    if sub is not None:
        ax.plot(sub['pct_cells_retained'], sub['edges_detected'] / 1e3, color='#555555',
                marker='s', ms=5, lw=1.5, zorder=2, label='random cell subsample\n(all 3 donors kept)')
    ax.scatter(lodo['pct_cells_retained'], lodo['edges_detected'] / 1e3, s=95,
               color=[COLORS[d] for d in lodo['left_out_donor']], zorder=3,
               edgecolor='white', linewidth=0.8, label='leave one donor out')
    for _, r in lodo.iterrows():
        ax.annotate(f"without {r['left_out_donor']}",
                    (r['pct_cells_retained'], r['edges_detected'] / 1e3),
                    textcoords='offset points', xytext=(7, -2), fontsize=7.5)
    ax.scatter([100], [EDGES['full'] / 1e3], s=95, color='#D62728', marker='D', zorder=3)
    ax.annotate('full (3 donors)', (100, EDGES['full'] / 1e3), textcoords='offset points',
                xytext=(-8, -13), fontsize=7.5, ha='right')
    ax.set_xlabel('% of cells retained')
    ax.set_ylabel('Edges detected (thousands)')
    ax.set_xlim(45, 108)
    ax.set_ylim(0, 88)
    ax.legend(fontsize=6.5, frameon=False, loc='upper left')
    ax.set_title('C  Discoveries scale with cells retained', fontsize=9, fontweight='bold',
                 loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    fig.tight_layout()
    for ext in ['pdf', 'png']:
        out = os.path.join(args.output_dir, f"donor_contribution.{ext}")
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Wrote {out}")

    logger.info("All done!")


if __name__ == "__main__":
    main()
