# %%
"""
Log fold changes of the TYK2 and DEXI cascade genes across every run, with standard errors.

Plots the log fold changes themselves, one line per gene across the runs, with +/- 1 standard
error from DESeq2, rather than encoding effect size as a marker area that cannot be read
quantitatively.

The standard errors are the point of the figure. Two mechanisms can remove a call at FDR 0.1:

    the estimate moves towards zero   -> a genuine change in the effect, e.g. a donor-specific program
    the standard error grows          -> loss of power, the estimate is unchanged

Plotting LFC +/- lfcSE separates them directly, which the significance flags alone cannot do. Each panel
also reports the mean |LFC| and the mean lfcSE relative to the full run, so the two components can be
compared as numbers.

Input is tyk2_dexi_profiles_all.tsv: the full profiles (all 8,001 genes) for the perturbations
TYK2, DEXI and SSB, from all seven runs. It is extracted once with awk, which is far faster than
parsing the ~27 GB of chunk CSVs with pandas. Run this from /g/stegle/schrod/code/TCell:

    OUT=Donor/output/tyk2_dexi_profiles_all.tsv
    printf "run\\tcontrast\\tvariable\\tbaseMean\\tlog_fc\\tlfcSE\\tp_value\\tadj_p_value\\n" > $OUT
    for spec in "full:pseudobulk_deseq2_chunks50_min3bulks_ribosomal" \\
                "A:Donor/pseudobulk_deseq2_without_donor_A" \\
                "B:Donor/pseudobulk_deseq2_without_donor_B" \\
                "C:Donor/pseudobulk_deseq2_without_donor_C" \\
                "keep90:Subsample/pseudobulk_deseq2_keep90_seed0" \\
                "keep75:Subsample/pseudobulk_deseq2_keep75_seed0" \\
                "keep50:Subsample/pseudobulk_deseq2_keep50_seed0"; do
      run="${spec%%:*}"; dir="${spec#*:}"
      awk -F, -v R="$run" 'FNR>1 && ($9=="TYK2" || $9=="DEXI" || $9=="SSB") \\
        {print R"\\t"$9"\\t"$2"\\t"$3"\\t"$4"\\t"$5"\\t"$7"\\t"$8}' "$dir"/*.csv >> $OUT
    done

The chunk CSV columns are index, variable, baseMean, log_fc, lfcSE, stat, p_value, adj_p_value,
contrast -- hence the $9/$2/$3/$4/$5/$7/$8 selection above.

Outputs (to --output-dir):
    tyk2_dexi_lfc_table.tsv            LFC, lfcSE and adj_p per gene and run
    tyk2_dexi_lfc_components.tsv       mean |LFC| and mean lfcSE per cascade step and run
    tyk2_dexi_lfc_across_runs.pdf/.png
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

FDR = 0.1
RUNS = ['full', 'A', 'B', 'C', 'keep90', 'keep75', 'keep50']
RUN_LABELS = {'full': 'full\n(3 donors)', 'A': 'without A\n(79%)', 'B': 'without B\n(65%)',
              'C': 'without C\n(56%)', 'keep90': 'random\n90%', 'keep75': 'random\n75%',
              'keep50': 'random\n50%'}

ISG = ['MX1', 'OAS1', 'OAS2', 'ISG15', 'IFI6', 'IFI44', 'IFI44L', 'IFIH1']
SSB_TARGETS = ['SLC7A5', 'SLC1A5', 'SLC3A2', 'WARS', 'HARS', 'EPRS']
MAX_DERIVED = 16

C_LODO_SPAN = '#FDE9E9'      # shading behind the leave-one-donor-out columns
C_MEAN = '#111111'


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor/output")
    return p.parse_args()


def lfc_panel(ax, d, contrast, genes, title, runs):
    """One line per gene: LFC across runs, with +/- 1 lfcSE and filled markers where significant."""
    sub = d[d.contrast == contrast]
    lfc = sub.pivot(index='variable', columns='run', values='log_fc')
    se = sub.pivot(index='variable', columns='run', values='lfcSE')
    ap = sub.pivot(index='variable', columns='run', values='adj_p_value')
    genes = [g for g in genes if g in lfc.index]

    x = np.arange(len(runs))
    # shade the leave-one-donor-out block so it reads apart from the random subsamples
    n_lodo = sum(1 for r in runs if r in ('A', 'B', 'C'))
    if n_lodo:
        ax.axvspan(0.5, 0.5 + n_lodo, color=C_LODO_SPAN, zorder=0)
    ax.axhline(0, color='#666666', lw=0.8, zorder=1)

    cmap = plt.get_cmap('tab20')
    for gi, g in enumerate(genes):
        y = np.array([lfc.loc[g, r] for r in runs], dtype=float)
        e = np.array([se.loc[g, r] for r in runs], dtype=float)
        col = cmap(gi % 20)
        ax.errorbar(x, y, yerr=e, color=col, lw=1.1, alpha=0.85, zorder=2,
                    elinewidth=0.8, capsize=2, marker='None')
        sig = np.array([ap.loc[g, r] < FDR for r in runs])
        ax.scatter(x[sig], y[sig], s=26, color=col, zorder=4, edgecolor='white', linewidth=0.5)
        ax.scatter(x[~sig], y[~sig], s=26, facecolor='white', edgecolor=col, linewidth=1.1, zorder=4)
        if len(genes) <= 8:
            ax.annotate(g, (x[0], y[0]), textcoords='offset points', xytext=(-6, 0),
                        ha='right', va='center', fontsize=6.5, color=col)

    # mean |LFC| across the gene set, as a heavy line
    m = np.array([lfc.loc[genes, r].abs().mean() for r in runs])
    sgn = np.sign(np.nanmean(lfc.loc[genes, 'full'].values))
    ax.plot(x, sgn * m, color=C_MEAN, lw=2.2, ls='--', zorder=5,
            label=f'mean |LFC| (signed)')

    ax.set_xticks(x)
    ax.set_xticklabels([RUN_LABELS[r] for r in runs], fontsize=6.5, rotation=35, ha='right')
    ax.set_xlim(-0.6, len(runs) - 0.4)
    ax.set_ylabel('log2 fold change')
    ax.set_title(title, fontsize=8.5, fontweight='bold', loc='left')
    ax.legend(fontsize=6.5, frameon=False, loc='lower right')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    return lfc, se, ap, genes


def main():
    args = parse_args()
    od = args.output_dir
    prof = os.path.join(od, "tyk2_dexi_profiles_all.tsv")
    if not os.path.exists(prof):
        logger.error(f"{prof} not found. Extract it with the awk command in the "
                     f"docstring of this file.")
        sys.exit(1)

    d = pd.read_csv(prof, sep='\t')
    runs = [r for r in RUNS if r in set(d['run'])]
    logger.info(f"Runs: {', '.join(runs)}; {len(d):,} rows")

    dexi_hits = (d[(d.contrast == 'DEXI') & (d.variable != 'DEXI') & (d.run == 'full') &
                   (d.adj_p_value < FDR)].sort_values('adj_p_value')['variable']
                 .head(MAX_DERIVED).tolist())

    steps = [('TYK2', ISG, 'A  TYK2 -> interferon program (Fig. 5B)'),
             ('DEXI', dexi_hits, 'B  DEXI -> its full-run targets (Fig. 5C, hop 1)'),
             ('SSB', SSB_TARGETS, 'C  SSB -> AA transport / tRNA (Fig. 5C, hop 2)')]

    # %%
    fig, axes = plt.subplots(1, 4, figsize=(18, 4.2),
                             gridspec_kw={'wspace': 0.3})
    comp_rows, tbl_rows = [], []
    for ax, (ctr, genes, title) in zip(axes, steps):
        lfc, se, ap, genes = lfc_panel(ax, d, ctr, genes, title, runs)
        for r in runs:
            comp_rows.append({
                'step': title.split('  ', 1)[1], 'contrast': ctr, 'run': r, 'n_genes': len(genes),
                'mean_abs_lfc': round(float(lfc.loc[genes, r].abs().mean()), 3),
                'mean_lfcSE': round(float(se.loc[genes, r].mean()), 3),
                'n_significant': int((ap.loc[genes, r] < FDR).sum()),
            })
            for g in genes:
                tbl_rows.append({'contrast': ctr, 'gene': g, 'run': r,
                                 'log_fc': lfc.loc[g, r], 'lfcSE': se.loc[g, r],
                                 'adj_p_value': ap.loc[g, r],
                                 'significant': bool(ap.loc[g, r] < FDR)})

    # %%
    # Panel D: the on-target knockdowns. A perturbation that was knocked down less well in one run
    # would produce a weaker downstream response for purely technical reasons, so this is the control
    # that has to be read alongside panels A-C.
    ax = axes[3]
    x = np.arange(len(runs))
    n_lodo = sum(1 for r in runs if r in ('A', 'B', 'C'))
    if n_lodo:
        ax.axvspan(0.5, 0.5 + n_lodo, color=C_LODO_SPAN, zorder=0)
    ax.axhline(0, color='#666666', lw=0.8, zorder=1)
    on_colors = {'TYK2': '#D62728', 'DEXI': '#2C7FB8', 'SSB': '#18974C'}
    for ctr, col in on_colors.items():
        s = d[(d.contrast == ctr) & (d.variable == ctr)].set_index('run')
        y = np.array([s.loc[r, 'log_fc'] for r in runs], dtype=float)
        e = np.array([s.loc[r, 'lfcSE'] for r in runs], dtype=float)
        sig = np.array([s.loc[r, 'adj_p_value'] < FDR for r in runs])
        bm = float(s['baseMean'].iloc[0])
        ax.errorbar(x, y, yerr=e, color=col, lw=1.4, zorder=2, elinewidth=0.9, capsize=2)
        ax.scatter(x[sig], y[sig], s=34, color=col, zorder=4, edgecolor='white', linewidth=0.5)
        ax.scatter(x[~sig], y[~sig], s=34, facecolor='white', edgecolor=col, linewidth=1.2, zorder=4)
        ax.annotate(f'{ctr} (baseMean {bm:.0f})', (x[-1], y[-1]), textcoords='offset points',
                    xytext=(6, 0), ha='left', va='center', fontsize=6.5, color=col)
        for r in runs:
            tbl_rows.append({'contrast': ctr, 'gene': f'{ctr} (on-target)', 'run': r,
                             'log_fc': s.loc[r, 'log_fc'], 'lfcSE': s.loc[r, 'lfcSE'],
                             'adj_p_value': s.loc[r, 'adj_p_value'],
                             'significant': bool(s.loc[r, 'adj_p_value'] < FDR)})
    ax.set_xticks(x)
    ax.set_xticklabels([RUN_LABELS[r] for r in runs], fontsize=6.5, rotation=35, ha='right')
    ax.set_xlim(-0.6, len(runs) + 0.6)
    ax.set_ylabel('log2 fold change')
    ax.set_title('D  On-target knockdown (control)', fontsize=8.5, fontweight='bold', loc='left')
    for s_ in ['top', 'right']:
        ax.spines[s_].set_visible(False)

    # TYK2's interferon response normalised by how well TYK2 itself was knocked down
    ot = d[(d.contrast == 'TYK2') & (d.variable == 'TYK2')].set_index('run')
    isg = d[(d.contrast == 'TYK2') & (d.variable.isin(ISG))]
    m_isg = isg.groupby('run')['log_fc'].apply(lambda v: v.abs().mean())
    norm = pd.DataFrame({
        'run': runs,
        'on_target_abs_lfc': [abs(ot.loc[r, 'log_fc']) for r in runs],
        'on_target_lfcSE': [ot.loc[r, 'lfcSE'] for r in runs],
        'isg_mean_abs_lfc': [m_isg[r] for r in runs]})
    norm['isg_per_unit_knockdown'] = (norm['isg_mean_abs_lfc'] / norm['on_target_abs_lfc']).round(3)
    base = float(norm.loc[norm.run == 'full', 'isg_per_unit_knockdown'].iloc[0])
    norm['ratio_vs_full'] = (norm['isg_per_unit_knockdown'] / base).round(2)
    norm.to_csv(os.path.join(od, "tyk2_isg_normalised_by_knockdown.tsv"), sep='\t', index=False)
    logger.info("\nTYK2 interferon response normalised by on-target knockdown:\n" +
                norm.round(3).to_string(index=False))
    logger.info("  Caution: the denominator is itself a noisy estimate (TYK2 baseMean ~23), so this "
                "ratio is only interpretable together with the raw values in panels A and D.")

    comp = pd.DataFrame(comp_rows)
    for c in ['mean_abs_lfc', 'mean_lfcSE']:
        base = comp.merge(comp[comp.run == 'full'][['contrast', c]], on='contrast',
                          suffixes=('', '_full'))[f'{c}_full']
        comp[f'{c}_vs_full'] = (comp[c].values / base.values).round(3)
    comp.to_csv(os.path.join(od, "tyk2_dexi_lfc_components.tsv"), sep='\t', index=False)
    pd.DataFrame(tbl_rows).to_csv(os.path.join(od, "tyk2_dexi_lfc_table.tsv"), sep='\t', index=False)

    for ctr in comp.contrast.unique():
        s = comp[comp.contrast == ctr]
        logger.info(f"\n{ctr}: mean |LFC| and mean lfcSE relative to the full run\n" +
                    s[['run', 'n_significant', 'mean_abs_lfc', 'mean_abs_lfc_vs_full',
                       'mean_lfcSE', 'mean_lfcSE_vs_full']].to_string(index=False))

    handles = [plt.Line2D([], [], marker='o', ls='', mfc='#555555', mec='#555555', ms=5,
                          label=f'significant (FDR < {FDR})'),
               plt.Line2D([], [], marker='o', ls='', mfc='white', mec='#555555', ms=5,
                          label='not significant'),
               plt.Line2D([], [], color='#555555', lw=1, label='error bars: +/- 1 lfcSE'),
               plt.Rectangle((0, 0), 1, 1, fc=C_LODO_SPAN, label='leave one donor out')]
    fig.legend(handles=handles, fontsize=7, frameon=False, ncol=4,
               loc='lower center', bbox_to_anchor=(0.5, -0.13))

    for ext in ['pdf', 'png']:
        out = os.path.join(od, f"tyk2_dexi_lfc_across_runs.{ext}")
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Wrote {out}")

    logger.info("All done!")


if __name__ == "__main__":
    main()
