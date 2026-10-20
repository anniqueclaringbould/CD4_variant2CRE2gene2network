"""Effect estimates and their significance across donor and random removal, as one heatmap.

Each column is one of the genes a cascade node significantly moves in the full three-donor
analysis; each row is one analysis. Colour carries the log2 fold change and a separate bar
carries significance, so a gene that keeps its effect but loses its FDR call -- the pattern
donor removal produces -- is visible as a row that keeps its colour and loses its bar.

Reads the per-gene tables written by 03_11_tyk2_dexi_lfc_across_runs.py.

Outputs (to --output-dir):
    node_effect_significance_heatmap.pdf / .png
"""

from __future__ import annotations

import argparse
import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import Patch, Rectangle

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 8

FDR = 0.1
NODES = ['DEXI', 'SSB', 'TYK2']
COLOR = {'DEXI': '#2C7FB8', 'SSB': '#18974C', 'TYK2': '#E67E22'}
RUNS = ['A', 'B', 'C', 'keep90', 'keep75', 'keep50']
LABEL = {'full': 'full', 'A': 'without A', 'B': 'without B', 'C': 'without C',
         'keep90': 'random 90%', 'keep75': 'random 75%', 'keep50': 'random 50%'}
RANDOM_REMOVED = {'keep90': 10.0, 'keep75': 25.0, 'keep50': 50.0}

#: The named programmes from the manuscript's locus figures, read from the curated table rather
#: than restated here. Within a node the programme genes are shown first, so a reader can see
#: whether the programme behaves like the rest of that node's downstream set or differently.
PROGRAMME_LABEL = {'DEXI': 'targets', 'SSB': 'AA transport / tRNA', 'TYK2': 'interferon'}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--output-dir', default="/g/stegle/schrod/code/TCell/Donor/output")
    return p.parse_args()


def save(fig, od, stem):
    for ext in ('pdf', 'png'):
        p = os.path.join(od, f"{stem}.{ext}")
        fig.savefig(p, dpi=300, bbox_inches='tight')
    print(f"Wrote {os.path.join(od, stem)}.pdf/.png")


def load(od):
    prof = pd.read_csv(os.path.join(od, "tyk2_dexi_profiles_all.tsv"), sep='\t')
    sets, frames = {}, []
    for node in NODES:
        d = prof[prof['contrast'] == node]
        full = d[d['run'] == 'full'].set_index('variable')
        tgt = full[(full['adj_p_value'] < FDR) & (full.index != node)].index
        sets[node] = tgt
        for run in ['full'] + RUNS:
            s = d[d['run'] == run].set_index('variable').reindex(tgt)
            frames.append(pd.DataFrame({
                'node': node, 'run': run, 'gene': tgt,
                'lfc': s['log_fc'].values, 'se': s['lfcSE'].values,
                'sig': (s['adj_p_value'] < FDR).values,
                'lfc_full': full.loc[tgt, 'log_fc'].values,
            }))
    return sets, pd.concat(frames, ignore_index=True)


def heatmap(sets, t, od, stem='node_effect_significance_heatmap'):
    """Genes x runs: colour is the log2 fold change, the bar beneath each row is significance.

    Within a node the genes of the named programme come first, marked in the node's colour in the
    strip above the panel, then the remaining downstream genes in grey.
    """
    prog_tbl = pd.read_csv(os.path.join(od, "tyk2_dexi_lfc_table.tsv"), sep='\t')
    prog = {n: set(prog_tbl[(prog_tbl['contrast'] == n)
                            & (~prog_tbl['gene'].str.contains('on-target'))]['gene'])
            for n in NODES}

    order, bounds, subs = [], [], []
    for node in NODES:
        f = t[(t['node'] == node) & (t['run'] == 'full')].copy()
        f['in_prog'] = ~f['gene'].isin(prog[node])        # programme genes first
        start = len(order)
        for is_other, label in ((False, PROGRAMME_LABEL[node]), (True, 'other')):
            g = f[f['in_prog'] == is_other].sort_values('lfc')['gene']
            if len(g):
                subs.append((label, len(order), len(order) + len(g), node))
                order += [(node, x) for x in g]
        bounds.append((node, start, len(order)))

    cols = ['full'] + RUNS
    M = np.full((len(order), len(cols)), np.nan)
    S = np.zeros_like(M, dtype=bool)
    idx = {(n, g): i for i, (n, g) in enumerate(order)}
    for r in t.itertuples():
        if r.run in cols:
            i, j = idx[(r.node, r.gene)], cols.index(r.run)
            M[i, j] = r.lfc
            S[i, j] = bool(r.sig)

    # Transposed: runs down, genes across. 158 genes is a long axis and 7 runs a short one, so
    # this fits a page width and leaves the run labels horizontal.
    M, S = M.T, S.T

    # Significance as a separate barcode under each run rather than as a wash over it. Fading the
    # colour made the two readings compete: a pale cell could be a small effect or a
    # non-significant one, and the magnitude -- the thing the figure is actually about -- lost.
    v = np.nanmax(np.abs(M))
    norm = TwoSlopeNorm(0, -v, v)
    rgba = plt.get_cmap('RdBu_r')(norm(M))
    rgba[np.isnan(M)] = (1, 1, 1, 1)

    BAND, BAR = 6, 1                      # rows of colour, then one thin barcode row
    H = BAND + BAR
    img = np.ones((len(cols) * H, M.shape[1], 4))
    for k in range(len(cols)):
        img[k * H:k * H + BAND] = rgba[k]
        bar = np.ones((M.shape[1], 4))
        bar[S[k]] = (0.15, 0.15, 0.15, 1)
        img[k * H + BAND] = bar

    # Two columns: heatmap and a dedicated axis for the colourbar. The per-run correlations with
    # the full analysis live in their own figure, so repeating them here only competes with the
    # two things this one is for -- the effect pattern and the significance barcode.
    fig, (ax, axc) = plt.subplots(
        1, 2, figsize=(7.4, 2.6),
        gridspec_kw={'width_ratios': [8, 0.18], 'wspace': 0.03})
    ax.imshow(img, aspect='auto', interpolation='nearest')
    ax.set_yticks([k * H + (BAND - 1) / 2 for k in range(len(cols))])
    ax.set_yticklabels([LABEL[c] for c in cols], fontsize=7)
    ax.set_xticks([])
    ax.tick_params(axis='y', length=0)
    for node, lo, hi in bounds:
        if hi < len(order):
            ax.axvline(hi - 0.5, color='black', lw=0.9)
        ax.text((lo + hi) / 2 - 0.5, -3.6, f"{node} ({hi - lo})", ha='center', va='bottom',
                fontsize=7.5, color=COLOR[node], fontweight='bold')
    # Programme membership as a strip above the panel rather than as inline labels: the SSB
    # programme is five columns wide and no horizontal label fits over it.
    for label, lo, hi, node in subs:
        is_prog = label != 'other'
        if lo > 0 and not any(lo == b[1] for b in bounds):
            ax.axvline(lo - 0.5, color='black', lw=0.6, ls=(0, (2, 1.6)))
        ax.add_patch(Rectangle((lo - 0.5, -3.1), hi - lo, 1.5,
                               facecolor=COLOR[node] if is_prog else '#d9d9d9',
                               edgecolor='none', clip_on=False))

    # The barcode is explained by a key above the panel rather than by a sentence beneath it,
    # so the figure carries no caption text of its own.
    ax.legend(handles=[Patch(facecolor='0.15', edgecolor='none',
                             label='significant at FDR 0.1'),
                       Patch(facecolor='#7a7a7a', edgecolor='none',
                             label='named programme (coloured) vs other downstream genes (grey)')],
              ncol=2, columnspacing=1.6,
              loc='upper center', bbox_to_anchor=(0.5, -0.03), frameon=False, fontsize=7,
              handlelength=1.1, handleheight=0.7, handletextpad=0.5, borderpad=0)

    sm = plt.cm.ScalarMappable(norm=norm, cmap='RdBu_r')
    cb = fig.colorbar(sm, cax=axc)
    cb.set_label('log$_2$FC', fontsize=7.5)
    cb.ax.tick_params(labelsize=7)
    fig.tight_layout()
    save(fig, od, stem)


def main() -> int:
    a = parse_args()
    sets, t = load(a.output_dir)
    print({n: len(s) for n, s in sets.items()})
    heatmap(sets, t, a.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
