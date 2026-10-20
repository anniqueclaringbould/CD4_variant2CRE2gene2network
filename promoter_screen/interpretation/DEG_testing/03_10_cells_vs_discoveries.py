# %%
"""
Cells contributed by each donor against the number of links discovered when that donor is left out.

The two are expected to be anticorrelated: the more a donor contributes, the more discoveries are lost
by dropping it. This is the compact statement of why the leave-one-donor-out (LODO) runs differ so
much, since the donors contribute 21%, 35% and 44% of the cells rather than a third each.

With three donors the correlation is descriptive only, not inferential, and is labelled as such on the
plot. n = 3 points cannot support a p value.

Reads the tables written by 03_07 (donor_contribution.tsv) and 03_05 (lodo_summary.tsv), so all
numbers come from the analyses rather than being restated here.

Outputs (to --output-dir):
    cells_vs_discoveries.tsv
    cells_vs_discoveries.pdf / .png
"""

import os
import sys
import argparse

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

COLORS = {'A': '#2C7FB8', 'B': '#18974C', 'C': '#E67E22'}
C_CELLS = '#7F7F7F'
C_EDGES = '#D62728'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor/output")
    return p.parse_args()


def main():
    args = parse_args()
    od = args.output_dir
    f_donor = os.path.join(od, "donor_contribution.tsv")
    f_lodo = os.path.join(od, "lodo_summary.tsv")
    for f in [f_donor, f_lodo]:
        if not os.path.exists(f):
            logger.error(f"{f} not found, run 03_07 and 03_05 first.")
            sys.exit(1)

    donor = pd.read_csv(f_donor, sep='\t', index_col=0)
    lodo = pd.read_csv(f_lodo, sep='\t')

    d = pd.DataFrame({
        'donor': donor.index,
        'cells_in_donor': donor['total_cells'].values,
        'pct_of_all_cells': donor['pct_of_all_cells'].values,
    }).merge(lodo[['left_out_donor', 'edges_detected', 'perturbations_tested']],
             left_on='donor', right_on='left_out_donor').drop(columns='left_out_donor')
    d = d.sort_values('cells_in_donor').reset_index(drop=True)
    d['pct_of_full_run_edges'] = (100 * d['edges_detected'] / lodo['edges_detected'].max()).round(1)
    d.to_csv(os.path.join(od, "cells_vs_discoveries.tsv"), sep='\t', index=False)

    r = float(np.corrcoef(d['cells_in_donor'], d['edges_detected'])[0, 1])
    logger.info("\n" + d.to_string(index=False))
    logger.info(f"\nPearson r (cells in donor vs edges when left out) = {r:.3f}  "
                f"[descriptive only, n={len(d)}]")

    # random cell subsampling reference from 03_08. Expressed as cells *removed*, so it shares the
    # x axis with the LODO points and answers whether removing a donor costs more than removing the
    # same number of random cells.
    sub_f = os.path.join(od, "subsample_summary.tsv")
    sub = None
    if os.path.exists(sub_f):
        total_cells = int(pd.read_csv(f_donor, sep='\t')['total_cells'].sum())
        sub = pd.read_csv(sub_f, sep='\t')
        sub['cells_removed'] = (total_cells - sub['cells_retained']).astype(int)
        sub = pd.concat([pd.DataFrame([{'cells_removed': 0,
                                        'edges_detected': lodo['edges_detected'].max()}]), sub],
                        ignore_index=True).sort_values('cells_removed')
        logger.info("\nRandom subsampling reference (cells removed vs edges):\n" +
                    sub[['cells_removed', 'edges_detected']].to_string(index=False))
        for row in d.itertuples():
            exp = np.interp(row.cells_in_donor, sub['cells_removed'], sub['edges_detected'])
            logger.info(f"  without {row.donor}: {row.edges_detected:,} edges vs {exp:,.0f} "
                        f"expected when removing {row.cells_in_donor:,} random cells "
                        f"({100 * (row.edges_detected / exp - 1):+.1f}%)")
    else:
        logger.warning(f"{sub_f} not found, plotting without the random subsampling reference")

    # %%
    fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))

    # A, mirrored bars: what a donor contributes vs what is left when it is dropped
    ax = axes[0]
    x = np.arange(len(d))
    w = 0.38
    ax.bar(x - w / 2, d['cells_in_donor'] / 1e3, width=w, color=C_CELLS,
           label='cells contributed by donor')
    ax.set_ylabel('Cells contributed (thousands)', color=C_CELLS)
    ax.tick_params(axis='y', labelcolor=C_CELLS)
    ax.set_ylim(0, d['cells_in_donor'].max() / 1e3 * 1.25)

    ax2 = ax.twinx()
    ax2.bar(x + w / 2, d['edges_detected'] / 1e3, width=w, color=C_EDGES,
            label='links found without that donor')
    ax2.set_ylabel('Links discovered when left out (thousands)', color=C_EDGES)
    ax2.tick_params(axis='y', labelcolor=C_EDGES)
    ax2.set_ylim(0, d['edges_detected'].max() / 1e3 * 1.25)

    ax.set_xticks(x)
    ax.set_xticklabels([f"donor {row.donor}\n({row.pct_of_all_cells:.0f}% of cells)"
                        for row in d.itertuples()], fontsize=7.5)
    ax.set_title('A  More cells contributed, fewer links left',
                 fontsize=9, fontweight='bold', loc='left')
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=7, frameon=False, loc='upper center')
    for s in ['top']:
        ax.spines[s].set_visible(False)
        ax2.spines[s].set_visible(False)

    # B, the anticorrelation itself, against the random removal reference
    ax = axes[1]
    if sub is not None:
        ax.plot(sub['cells_removed'] / 1e3, sub['edges_detected'] / 1e3, color='#555555',
                marker='s', ms=5, lw=1.5, zorder=2,
                label='same number of random cells removed\n(all 3 donors kept)')
    ax.scatter(d['cells_in_donor'] / 1e3, d['edges_detected'] / 1e3, s=140,
               color=[COLORS[v] for v in d['donor']], zorder=3, edgecolor='white', linewidth=1,
               label='whole donor removed')
    for row in d.itertuples():
        ax.annotate(f"without {row.donor}",
                    (row.cells_in_donor / 1e3, row.edges_detected / 1e3),
                    textcoords='offset points', xytext=(0, 13), ha='center', fontsize=8)
    ax.set_xlabel('Cells removed (thousands)')
    ax.set_ylabel('Links discovered (thousands)')
    ax.set_ylim(0, max(d['edges_detected'].max(),
                       sub['edges_detected'].max() if sub is not None else 0) / 1e3 * 1.32)
    ax.legend(fontsize=6.5, frameon=False, loc='lower left')
    ax.set_title(f'B  Donor removal vs random removal\n'
                 f'(cells vs links when left out: r = {r:.2f}, n=3)',
                 fontsize=9, fontweight='bold', loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    fig.tight_layout()
    for ext in ['pdf', 'png']:
        out = os.path.join(od, f"cells_vs_discoveries.{ext}")
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Wrote {out}")

    logger.info("All done!")


if __name__ == "__main__":
    main()
