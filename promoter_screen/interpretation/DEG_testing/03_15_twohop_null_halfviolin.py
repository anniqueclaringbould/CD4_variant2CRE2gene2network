"""The two-hop matched-null result: the real network as a bar, the null as a half violin.

Three measures with different units -- a percentage of the transcriptome and two term counts --
share one axis, so every bar is annotated and the axis names both units. The null is drawn as a
density rather than a median and a whisker because the three nulls have quite different shapes,
and the hop-1 null in particular is a spike at zero with a thin tail that a summary statistic
hides. The sample is reflected about zero before the kernel density: without that the kernel
leaks mass below zero, which is impossible for a count, and the leak comes out of exactly the
boundary spike that matters.

Reads the per-draw table written by 03_14_twohop_expression_matched_null.py (twohop_null.tsv).

Outputs (to --output-dir):
    twohop_null_halfviolin.pdf / .png
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 8

RUN_DIR = "/g/stegle/schrod/code/TCell/Network/final_expr_1k"
#: Exact; rounding pct_hop2 to 56.0 before comparing moves the empirical p to 0.652.
REAL = {'pct_hop2': 56.047349, 'n_sig': 45, 'n_sig_hop2': 25}
PANEL = [
    ('pct_hop2', 'genes reached\nwithin 2 hops (%)', 'reach'),
    ('n_sig', 'enriched GO BP terms,\nhop 1 (815 genes)', 'hop 1'),
    ('n_sig_hop2', 'enriched GO BP terms,\nhop 2 (4,356 genes)', 'hop 2'),
]
#: Real in red, null in blue -- the same two colours used for the observed value and the
#: matched baseline in the other response figures.
REAL_C, NULL_C = '#D62728', '#2C7FB8'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--run-dir', default=RUN_DIR)
    p.add_argument('--output-dir', default=RUN_DIR)
    return p.parse_args()


def main() -> int:
    a = parse_args()
    t = pd.read_csv(os.path.join(a.run_dir, "twohop_null.tsv"), sep='\t')
    n = len(t)

    # Smooth density, with the sample reflected about zero before the KDE. Without that the
    # kernel leaks mass below zero -- impossible for a count -- and the leak is subtracted from
    # the spike at the boundary, which is exactly the feature of these nulls that matters.
    labels = ['genes reached\nwithin 2 hops (%)', 'enriched GO terms\nhop 1',
              'enriched GO terms\nhop 2']
    h2 = [plt.Rectangle((0, 0), 1, 1, color=REAL_C),
          plt.Rectangle((0, 0), 1, 1, color=NULL_C)]

    fig4, ax = plt.subplots(figsize=(5.8, 2.7))
    BAR_W, VIO_W = 0.40, 0.40
    ymax = 1.22 * max(max(np.percentile(t[c], 99.5), REAL[c]) for c, _, _ in PANEL)
    for i, (col, _, _) in enumerate(PANEL):
        real, med = REAL[col], t[col].median()
        ge = int((t[col] >= real).sum())
        pv = (ge + 1) / (n + 1)
        v = t[col].values.astype(float)

        # The bar's right edge and the violin's baseline are the same x, so the two touch.
        edge = float(i)
        bar_c = edge - BAR_W / 2
        ax.bar(bar_c, real, width=BAR_W, color=REAL_C, zorder=2)
        ax.text(bar_c, real + ymax * 0.012, f"{real:.0f}", ha='center', va='bottom',
                fontsize=8, color=REAL_C, fontweight='bold')

        grid = np.linspace(0, min(v.max(), ymax * 0.985), 512)
        dens = gaussian_kde(np.concatenate([v, -v]))(grid) * 2
        width = VIO_W * dens / dens.max()
        x0 = edge                           # flush with the bar's right edge
        ax.fill_betweenx(grid, x0, x0 + width, color=NULL_C, linewidth=0, zorder=2)
        w_med = float(np.interp(med, grid, width))
        ax.plot([x0, x0 + w_med], [med, med], color='#17455e', lw=1.4, zorder=3)
        ax.text(x0 + w_med + 0.03, med, f"{med:.0f}", ha='left', va='center',
                fontsize=8, color='#17455e')

        # p at a common height, so no label can collide with a tall density spire.
        ax.text(i, ymax * 0.93, f"$p$ = {pv:.3f}", ha='center', va='center', fontsize=7.5)

    ax.set_xticks(np.arange(len(PANEL)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel('per cent  /  no. of terms', fontsize=8)
    ax.set_ylim(0, ymax)
    ax.set_xlim(-0.62, len(PANEL) - 0.22)
    ax.tick_params(axis='y', labelsize=8)
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.legend(h2, ['CRE target network',
                   f'matched random networks, density and median (n={n:,})'],
              frameon=False, fontsize=7, loc='upper center', bbox_to_anchor=(0.5, -0.24),
              ncol=2, handletextpad=0.6)
    fig4.tight_layout()
    for ext in ('pdf', 'png'):
        out = os.path.join(a.output_dir, f"twohop_null_halfviolin.{ext}")
        fig4.savefig(out, dpi=300, bbox_inches='tight', pad_inches=0.08)
        print(f"Wrote {out}")

    for col, _, title in PANEL:
        ge = int((t[col] >= REAL[col]).sum())
        lo, hi = np.percentile(t[col], [5, 95])
        print(f"{title:6} real {REAL[col]:>5}  null median {t[col].median():>5.1f}  "
              f"5-95% {lo:.0f}-{hi:.0f}  {ge}/{n}  p={(ge + 1) / (n + 1):.3f}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
