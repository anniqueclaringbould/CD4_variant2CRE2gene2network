# %%
"""
Distribution of the cell state AUC scores, per cell and per pseudobulk.

Averaging cells into pseudobulks does two things to a score, and they pull in opposite directions:

    it removes within-bulk noise, so the between-bulk spread is much narrower than the between-cell
    spread, and the correlations between scores rise sharply

    it leaves a sampling error of roughly sd_cell / sqrt(n_cells), which for the many small pseudobulks
    is not negligible, so part of the between-bulk spread of a score is measurement error rather than
    biology

The second point matters for using these as covariates: a covariate measured with error is attenuated
towards zero (regression dilution), and adjusts for less than it appears to. This script quantifies both
effects and checks how much of each score is explained by donor, since a score that mostly tracks donor
is largely redundant with the donor term already in the design.

Outputs (to --output-dir):
    cell_score_distributions.tsv           per score: cell and pseudobulk moments, noise fraction, donor R2
    cell_score_distributions.pdf/.png      paired violins, noise vs bulk size, and the donor split
    cell_score_overlay.pdf/.png            per score, the cell and pseudobulk densities overlaid
    cell_score_donor_overlay.pdf/.png      per score, the pseudobulk distributions of the three donors
                                           overlaid
    cell_score_within_donor_corr.tsv       score correlations overall and within each donor
"""

import os
import sys
import logging
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

# the logger is defined here rather than taken from utils, because `from utils import *` pulls in
# scanpy -> leidenalg -> igraph, which this script does not need and which intermittently hangs on
# import. Same format as utils.get_logger so the log output is unchanged.
logging.basicConfig(format="[%(asctime)s] %(levelname)s:%(name)s: %(message)s", stream=sys.stdout)
logger = logging.getLogger('utils')
logger.setLevel(logging.INFO)

SCORE_FILE = ("/g/steinmetz/project/otar_2063/promoter_screen/interpretation/cell_scoring/"
              "output/Representative_AUC_scores_per_cell.tsv")
CELL_H5AD = "/g/stegle/schrod/data/T_Cell/promoter_single_guide_hvg_raw_counts.h5ad"
BULK_SCORES = "/g/stegle/schrod/data/T_Cell/pseudobulk_cell_scores.tsv"

C_CELL, C_BULK = '#B0B0B0', '#2C7FB8'
DONOR_COLORS = {'A': '#2C7FB8', 'B': '#18974C', 'C': '#E67E22'}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--score-file', type=str, default=SCORE_FILE)
    p.add_argument('--cell-h5ad', type=str, default=CELL_H5AD)
    p.add_argument('--bulk-scores', type=str, default=BULK_SCORES)
    p.add_argument('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor/output")
    return p.parse_args()


def read_obs(path, cols):
    out = {}
    with h5py.File(path, 'r') as f:
        idx_key = f['obs'].attrs.get('_index', '_index')
        index = pd.Index([x.decode() for x in f['obs'][idx_key][:]])
        for c in cols:
            g = f['obs'][c]
            if isinstance(g, h5py.Group):
                cats = g['categories'][:].astype(str)
                out[c] = pd.Categorical.from_codes(g['codes'][:], cats).astype(str)
            else:
                out[c] = g[:]
    return pd.DataFrame(out, index=index)


def main():
    args = parse_args()
    od = args.output_dir
    os.makedirs(od, exist_ok=True)

    # %%
    logger.info("Loading per-cell scores and cell metadata...")
    scores = pd.read_csv(args.score_file, sep='\t').set_index('cell')
    obs = read_obs(args.cell_h5ad, ['donor', 'guide'])
    cells = obs.join(scores, how='inner')
    raw_names = list(scores.columns)
    safe = {c: c.replace('.', '_') for c in raw_names}
    cells = cells.rename(columns=safe)
    terms = [safe[c] for c in raw_names]
    logger.info(f"  {len(cells):,} cells used by the pseudobulks, {len(terms)} scores")

    bulks = pd.read_csv(args.bulk_scores, sep='\t', index_col=0)
    bulks['donor'] = pd.Index(bulks.index).str.split('_').str[0]
    logger.info(f"  {len(bulks):,} pseudobulks")

    # %%
    # per score: cell and pseudobulk moments, the sampling noise floor, and the donor share
    rows = []
    for t in terms:
        sd_cell = float(cells[t].std())
        sd_bulk = float(bulks[t].std())
        # expected sd of a bulk mean from sampling alone, averaged over the observed bulk sizes
        exp_noise_var = float((sd_cell ** 2 / bulks['n_cells_scored']).mean())
        noise_frac = exp_noise_var / (sd_bulk ** 2)
        # share of the between-bulk variance explained by donor (one-way ANOVA R2)
        grp = bulks.groupby('donor')[t]
        ss_between = float(((grp.mean() - bulks[t].mean()) ** 2 * grp.size()).sum())
        ss_total = float(((bulks[t] - bulks[t].mean()) ** 2).sum())
        rows.append({
            'score': t,
            'cell_mean': round(float(cells[t].mean()), 3), 'cell_sd': round(sd_cell, 3),
            'cell_min': round(float(cells[t].min()), 3), 'cell_max': round(float(cells[t].max()), 3),
            'bulk_mean': round(float(bulks[t].mean()), 3), 'bulk_sd': round(sd_bulk, 3),
            'sd_ratio_bulk_over_cell': round(sd_bulk / sd_cell, 3),
            'sampling_noise_frac_of_bulk_var': round(noise_frac, 3),
            'donor_R2_of_bulk_var': round(ss_between / ss_total, 3),
        })
    tbl = pd.DataFrame(rows)
    tbl.to_csv(os.path.join(od, "cell_score_distributions.tsv"), sep='\t', index=False)
    logger.info("\n" + tbl.to_string(index=False))

    # %%
    fig = plt.figure(figsize=(13.5, 8))
    gs = fig.add_gridspec(2, 2, hspace=0.55, wspace=0.25, height_ratios=[1, 1])

    # A, cell vs pseudobulk distributions, paired
    ax = fig.add_subplot(gs[0, :])
    pos = np.arange(len(terms))
    vp1 = ax.violinplot([cells[t].values for t in terms], positions=pos - 0.19, widths=0.34,
                        showextrema=False, showmedians=True)
    vp2 = ax.violinplot([bulks[t].values for t in terms], positions=pos + 0.19, widths=0.34,
                        showextrema=False, showmedians=True)
    for b in vp1['bodies']:
        b.set_facecolor(C_CELL); b.set_alpha(0.85)
    for b in vp2['bodies']:
        b.set_facecolor(C_BULK); b.set_alpha(0.85)
    for vp in (vp1, vp2):
        vp['cmedians'].set_color('#111111'); vp['cmedians'].set_linewidth(1)
    ax.set_xticks(pos)
    ax.set_xticklabels(terms, rotation=20, ha='right')
    ax.set_ylabel('AUC score')
    ax.set_ylim(0, 1)
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=C_CELL, label=f'per cell (n={len(cells):,})'),
                       plt.Rectangle((0, 0), 1, 1, fc=C_BULK, label=f'per pseudobulk (n={len(bulks):,})')],
              fontsize=7.5, frameon=False, loc='upper right', ncol=2)
    ax.set_title('A  Averaging into pseudobulks collapses the spread',
                 fontsize=9, fontweight='bold', loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    # B, spread of the pseudobulk mean against the number of cells in the bulk
    ax = fig.add_subplot(gs[1, 0])
    bins = [1, 3, 5, 10, 20, 40, 80, 160, 10 ** 4]
    lab = ['1-2', '3-4', '5-9', '10-19', '20-39', '40-79', '80-159', '160+']
    grp = pd.cut(bulks['n_cells_scored'], bins=bins, right=False, labels=lab)
    cmap = plt.get_cmap('tab10')
    for i, t in enumerate(terms):
        sd_by_bin = bulks.groupby(grp, observed=True)[t].std()
        ax.plot(range(len(sd_by_bin)), sd_by_bin.values, marker='o', ms=3.5, lw=1.3,
                color=cmap(i % 10), label=t)
        if i == 0:
            # sampling-noise prediction for this score: sd_cell / sqrt(n)
            ncent = np.array([np.sqrt(bulks.loc[grp == l, 'n_cells_scored'].mean()) for l in sd_by_bin.index])
            ax.plot(range(len(sd_by_bin)), float(cells[t].std()) / ncent, ls='--', lw=1.4,
                    color='#111111', label=f'sd_cell/sqrt(n) for {t}')
    ax.set_xticks(range(len(lab)))
    ax.set_xticklabels(lab, rotation=30, ha='right')
    ax.set_xlabel('cells per pseudobulk')
    ax.set_ylabel('sd of the pseudobulk score')
    ax.legend(fontsize=6, frameon=False, ncol=2)
    ax.set_title('B  Small pseudobulks have noisier scores', fontsize=9, fontweight='bold', loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    # C, pseudobulk scores by donor
    ax = fig.add_subplot(gs[1, 1])
    donors = sorted(bulks['donor'].unique())
    w = 0.8 / len(donors)
    for di, dn in enumerate(donors):
        sub = bulks[bulks.donor == dn]
        bp = ax.boxplot([sub[t].values for t in terms],
                        positions=pos + (di - (len(donors) - 1) / 2) * w, widths=w * 0.85,
                        showfliers=False, patch_artist=True,
                        medianprops=dict(color='#111111', lw=0.8))
        for b in bp['boxes']:
            b.set_facecolor(DONOR_COLORS.get(dn, '#999999')); b.set_alpha(0.75)
    ax.set_xticks(pos)
    ax.set_xticklabels(terms, rotation=25, ha='right', fontsize=7)
    ax.set_ylabel('pseudobulk score')
    ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, fc=DONOR_COLORS.get(d, '#999999'),
                                     label=f'donor {d}') for d in donors],
              fontsize=7, frameon=False, ncol=3, loc='upper right')
    ax.set_title('C  Donor explains part of every score', fontsize=9, fontweight='bold', loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    for ext in ['pdf', 'png']:
        out = os.path.join(od, f"cell_score_distributions.{ext}")
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Wrote {out}")

    # %%
    # Second figure: the two distributions overlaid per score. Densities are normalised so that the
    # 1.4 M cells and the 61 k pseudobulks are on the same footing, and shared bins are used so the
    # two histograms are directly comparable.
    ncol = 4
    nrow = int(np.ceil(len(terms) / ncol))
    fig2, axes2 = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.7 * nrow))
    axes2 = np.atleast_1d(axes2).ravel()

    for i, t in enumerate(terms):
        ax = axes2[i]
        cv, bv = cells[t].values, bulks[t].values
        lo = float(min(cv.min(), bv.min()))
        hi = float(max(np.percentile(cv, 99.9), bv.max()))
        edges = np.linspace(lo, hi, 81)
        ctr = 0.5 * (edges[:-1] + edges[1:])

        hc, _ = np.histogram(cv, bins=edges, density=True)
        hb, _ = np.histogram(bv, bins=edges, density=True)
        ax.fill_between(ctr, hc, step='mid', color=C_CELL, alpha=0.75, lw=0, label='per cell')
        ax.fill_between(ctr, hb, step='mid', color=C_BULK, alpha=0.55, lw=0, label='per pseudobulk')
        ax.plot(ctr, hc, drawstyle='steps-mid', color='#7A7A7A', lw=0.9)
        ax.plot(ctr, hb, drawstyle='steps-mid', color=C_BULK, lw=1.1)
        ax.axvline(cv.mean(), color='#7A7A7A', ls='--', lw=1)
        ax.axvline(bv.mean(), color=C_BULK, ls='--', lw=1)

        row = tbl[tbl.score == t].iloc[0]
        ax.set_title(f"{t}\nsd {row.cell_sd:.3f} -> {row.bulk_sd:.3f} "
                     f"({row.sd_ratio_bulk_over_cell:.2f}x)", fontsize=8, fontweight='bold', loc='left')
        ax.set_xlabel('AUC score')
        ax.set_ylabel('density' if i % ncol == 0 else '')
        ax.set_xlim(lo, hi)
        for s in ['top', 'right']:
            ax.spines[s].set_visible(False)
        if i == 0:
            ax.legend(fontsize=6.5, frameon=False, loc='upper right')

    for j in range(len(terms), len(axes2)):
        axes2[j].set_visible(False)

    fig2.suptitle('Cell-level and pseudobulk-level score distributions, overlaid '
                  '(dashed lines = means, densities normalised)',
                  fontsize=9.5, fontweight='bold', x=0.01, ha='left')
    fig2.tight_layout(rect=[0, 0, 1, 0.95])
    for ext in ['pdf', 'png']:
        out = os.path.join(od, f"cell_score_overlay.{ext}")
        fig2.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Wrote {out}")

    # %%
    # Third figure: the same overlay idea as above, but stratified by donor. One panel per score, the
    # three donors overlaid, so a score whose spread is mostly between donors is obvious as three
    # displaced humps rather than one shared distribution.
    donors = sorted(bulks['donor'].unique())

    fig3, axes3 = plt.subplots(nrow, ncol, figsize=(3.2 * ncol, 2.7 * nrow))
    axes3 = np.atleast_1d(axes3).ravel()

    for i, t in enumerate(terms):
        ax = axes3[i]
        vals = {dn: bulks.loc[bulks.donor == dn, t].values for dn in donors}
        allv = bulks[t].values
        cellv = cells[t].values
        # range covers both levels, the cell distribution is the wider of the two
        lo = float(min(allv.min(), cellv.min()))
        hi = float(max(np.percentile(cellv, 99.5), allv.max()))
        edges = np.linspace(lo, hi, 61)
        ctr = 0.5 * (edges[:-1] + edges[1:])

        # all single cells pooled, as a grey backdrop for the donor curves
        hc, _ = np.histogram(cellv, bins=edges, density=True)
        ax.fill_between(ctr, hc, step='mid', color=C_CELL, alpha=0.55, lw=0,
                        label='all cells', zorder=0)
        ax.plot(ctr, hc, drawstyle='steps-mid', color='#8A8A8A', lw=0.9, zorder=1)

        for dn in donors:
            h, _ = np.histogram(vals[dn], bins=edges, density=True)
            col = DONOR_COLORS.get(dn, '#999999')
            ax.fill_between(ctr, h, step='mid', color=col, alpha=0.4, lw=0, label=f'donor {dn}',
                            zorder=2)
            ax.plot(ctr, h, drawstyle='steps-mid', color=col, lw=1.2, zorder=3)
            ax.axvline(vals[dn].mean(), color=col, ls='--', lw=1, zorder=4)

        # how much of the between-bulk variance is between donors (one-way ANOVA R2)
        grp = bulks.groupby('donor')[t]
        ss_b = float(((grp.mean() - bulks[t].mean()) ** 2 * grp.size()).sum())
        ss_t = float(((bulks[t] - bulks[t].mean()) ** 2).sum())
        ax.set_title(f"{t}\ndonor R2 = {ss_b / ss_t:.2f}", fontsize=8, fontweight='bold', loc='left')
        ax.set_xlabel('pseudobulk AUC score')
        ax.set_ylabel('density' if i % ncol == 0 else '')
        ax.set_xlim(lo, hi)
        for s in ['top', 'right']:
            ax.spines[s].set_visible(False)
        if i == 0:
            ax.legend(fontsize=6.5, frameon=False, loc='upper right')

    for j in range(len(terms), len(axes3)):
        axes3[j].set_visible(False)

    fig3.suptitle('Pseudobulk score distributions by donor, with all single cells pooled in the '
                  'background (dashed lines = donor means, densities normalised)',
                  fontsize=9.5, fontweight='bold', x=0.01, ha='left')
    fig3.tight_layout(rect=[0, 0, 1, 0.95])

    # keep the within-donor correlation table, it is the number that matters for the design and it
    # is not a plot
    corr_rows = []
    for lab, sub in [('pooled', bulks)] + [(f'donor {d}', bulks[bulks.donor == d]) for d in donors]:
        C = sub[terms].corr()
        for a in terms:
            for b in terms:
                if a < b:
                    corr_rows.append({'group': lab, 'score_a': a, 'score_b': b,
                                      'pearson_r': round(float(C.loc[a, b]), 3)})
    pd.DataFrame(corr_rows).to_csv(os.path.join(od, "cell_score_within_donor_corr.tsv"),
                                   sep='\t', index=False)

    cw = pd.DataFrame(corr_rows).pivot_table(index=['score_a', 'score_b'], columns='group',
                                             values='pearson_r')
    within = cw[[c for c in cw.columns if c.startswith('donor')]].mean(axis=1)
    logger.info("\nScore pair correlations, pooled vs the mean within donor:")
    cmp = pd.DataFrame({'pooled': cw['pooled'], 'mean_within_donor': within.round(3)})
    cmp['drop'] = (cmp['pooled'] - cmp['mean_within_donor']).round(3)
    logger.info("\n" + cmp.sort_values('pooled', ascending=False).head(12).to_string())
    logger.info(f"\nmean |r| off-diagonal: pooled {cw['pooled'].abs().mean():.3f}, "
                f"within donor {within.abs().mean():.3f}")

    for ext in ['pdf', 'png']:
        out = os.path.join(od, f"cell_score_donor_overlay.{ext}")
        fig3.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Wrote {out}")

    logger.info("All done!")


if __name__ == "__main__":
    main()
