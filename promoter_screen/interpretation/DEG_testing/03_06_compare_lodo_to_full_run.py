# %%
"""
Compare the leave-one-donor-out (LODO) DESeq2 runs against the full three donor run, and quantify
log fold change agreement rather than only set overlap.

Two questions:
  1. Are the LODO edge sets subsets of the published full run, and what fraction of the published
     edges does each LODO run recover?
  2. Do the log fold changes agree across runs? Set overlap is harsh on effects sitting near the
     FDR threshold, so LFC correlation is the fairer measure of agreement.

As in 03_05, everything is restricted to perturbations tested in *all* runs, so that differences in
testing coverage (targets dropped for too few pseudobulks) cannot show up as disagreement.

To get unbiased LFCs the script does a second streaming pass over all four result directories,
extracting every requested (perturbation, gene) pair regardless of its p value. Three pair sets are
evaluated separately, because the selection determines what the correlation means:
  full_run_significant  edges called in the full run    "do the LODO runs reproduce published effects"
  union_detected        edges called in any run         selected on significance, so upward biased
  random               random commonly tested pairs     unselected baseline, mostly null

Outputs (to --output-dir):
    lodo_vs_full_lfc_matrix.parquet   per pair LFC and adj_p in all four runs (intermediate)
    lodo_vs_full_recovery.tsv         recovery and subset statistics per LODO condition
    lodo_vs_full_correlations.tsv     Pearson and Spearman per run pair and pair set
    lodo_vs_full_comparison.pdf/.png  recovery, correlation heatmap and LFC scatters
"""

import os
import sys
import argparse
import itertools
from functools import reduce

import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Illustrator-compatible, as in 04_02_Figure4E.py
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.size'] = 8

os.chdir("/g/stegle/schrod/code/TCell")
from utils import *

logger = get_logger(__name__)

FDR = 0.1
FULL = 'full'            # label for the three donor run
N_RANDOM = 200_000       # random commonly tested pairs for the unselected baseline
RANDOM_SEED = 42

C_BAR = '#2C7FB8'
C_ALT = '#18974C'
C_GREY = '#B0B0B0'


def scan_dir(deg_dir, label, cache_p=0.5):
    """Pass 1: tested contrasts, tested genes and near significant edges of one run."""
    files = sorted(f for f in os.listdir(deg_dir) if f.endswith('.csv'))
    logger.info(f"  {label}: reading {len(files)} chunk files")
    edges, contrasts, genes = [], set(), set()
    for i, fn in enumerate(files, start=1):
        df = pd.read_csv(os.path.join(deg_dir, fn),
                         usecols=['contrast', 'variable', 'log_fc', 'adj_p_value'])
        contrasts.update(df['contrast'].unique())
        genes.update(df['variable'].unique())
        edges.append(df[df['adj_p_value'] < cache_p])
        if i % 50 == 0 or i == len(files):
            logger.info(f"    {i}/{len(files)} files")
    edges = pd.concat(edges, ignore_index=True)
    edges['run'] = label
    return edges, contrasts, genes


def extract_pairs(deg_dir, label, target_pairs):
    """Pass 2: LFC and adj_p for a fixed set of pairs, regardless of significance."""
    files = sorted(f for f in os.listdir(deg_dir) if f.endswith('.csv'))
    logger.info(f"  {label}: extracting {len(target_pairs):,} pairs from {len(files)} files")
    out = []
    for i, fn in enumerate(files, start=1):
        df = pd.read_csv(os.path.join(deg_dir, fn),
                         usecols=['contrast', 'variable', 'log_fc', 'adj_p_value'])
        df['edge'] = df['contrast'] + '->' + df['variable']
        out.append(df.loc[df['edge'].isin(target_pairs), ['edge', 'log_fc', 'adj_p_value']])
        if i % 50 == 0 or i == len(files):
            logger.info(f"    {i}/{len(files)} files, {sum(len(o) for o in out):,} rows kept")
    res = pd.concat(out, ignore_index=True)
    if res['edge'].duplicated().any():
        logger.error(f"{label}: duplicated (perturbation, gene) pairs found, cannot align runs.")
        sys.exit(1)
    return res.rename(columns={'log_fc': f'lfc_{label}', 'adj_p_value': f'adjp_{label}'}) \
              .set_index('edge')


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lodo-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor",
                   help='Directory holding the pseudobulk_deseq2_without_donor_* directories')
    p.add_argument('--full-dir', type=str,
                   default="/g/stegle/schrod/code/TCell/pseudobulk_deseq2_chunks50_min3bulks_ribosomal",
                   help='Directory with the published three donor DESeq2 results')
    p.add_argument('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor/output",
                   help='Directory for the summary tables and figures')
    p.add_argument('--force-reload', action='store_true', help='Ignore cached intermediates')
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    lfc_file = os.path.join(args.output_dir, "lodo_vs_full_lfc_matrix.parquet")

    prefix = "pseudobulk_deseq2_without_donor_"
    lodo = sorted(d.replace(prefix, "") for d in os.listdir(args.lodo_dir) if d.startswith(prefix))
    runs = [FULL] + lodo
    dirs = {FULL: args.full_dir, **{c: os.path.join(args.lodo_dir, prefix + c) for c in lodo}}
    labels = {FULL: 'full (3 donors)', **{c: f'without donor {c}' for c in lodo}}
    logger.info(f"Runs: {', '.join(labels[r] for r in runs)}")

    # %%
    if os.path.exists(lfc_file) and not args.force_reload:
        logger.info(f"Loading cached LFC matrix from {lfc_file}")
        wide = pd.read_parquet(lfc_file)
    else:
        # Pass 1: significant edges and testing coverage of every run. The LODO runs were already
        # scanned by 03_05, so reuse those caches and only scan the full run when possible.
        logger.info("Pass 1: scanning runs for detected edges and testing coverage...")
        near, tested_contrasts, gene_sets = {}, {}, {}

        lodo_cache = [os.path.join(args.output_dir, f) for f in
                      ["lodo_edges_near_significant.parquet", "lodo_tested_contrasts.csv",
                       "lodo_tested_genes.csv"]]
        if all(os.path.exists(f) for f in lodo_cache):
            logger.info("  reusing the LODO caches written by 03_05")
            cached_edges = pd.read_parquet(lodo_cache[0]).rename(columns={'condition': 'run'})
            cached_contrasts = pd.read_csv(lodo_cache[1])
            cached_genes = pd.read_csv(lodo_cache[2])
            for c in lodo:
                near[c] = cached_edges[cached_edges['run'] == c]
                tested_contrasts[c] = set(cached_contrasts.loc[cached_contrasts[c], 'contrast'])
                gene_sets[c] = set(cached_genes.loc[cached_genes[c], 'gene'])
            to_scan = [FULL]
        else:
            to_scan = runs

        for r in to_scan:
            near[r], tested_contrasts[r], gene_sets[r] = scan_dir(dirs[r], r)

        common_contrasts = reduce(set.intersection, tested_contrasts.values())
        common_genes = reduce(set.intersection, gene_sets.values())
        logger.info("Tested perturbations: " +
                    ", ".join(f"{labels[r]}={len(tested_contrasts[r]):,}" for r in runs))
        logger.info(f"Tested in all runs: {len(common_contrasts):,} perturbations, "
                    f"{len(common_genes):,} genes")

        # detected edges per run, on common ground, self edges excluded
        sig_sets = {}
        for r in runs:
            d = near[r]
            d = d[d['contrast'].isin(common_contrasts) & d['variable'].isin(common_genes)]
            d = d[(d['adj_p_value'] < FDR) & (d['contrast'] != d['variable'])]
            sig_sets[r] = set(d['contrast'] + '->' + d['variable'])
            logger.info(f"  {labels[r]}: {len(sig_sets[r]):,} edges at FDR < {FDR}")

        # %%
        # Pair sets to evaluate
        union_detected = reduce(set.union, sig_sets.values())
        rng = np.random.default_rng(RANDOM_SEED)
        cl, gl = sorted(common_contrasts), sorted(common_genes)
        rand_pairs = {f"{cl[i]}->{gl[j]}" for i, j in
                      zip(rng.integers(0, len(cl), N_RANDOM), rng.integers(0, len(gl), N_RANDOM))}
        rand_pairs = {p for p in rand_pairs if p.split('->')[0] != p.split('->')[1]}
        target_pairs = union_detected | rand_pairs
        logger.info(f"Pair sets: union_detected={len(union_detected):,}, "
                    f"random={len(rand_pairs):,}, total to extract={len(target_pairs):,}")

        # %%
        # Pass 2: LFC and adj_p for those pairs in every run
        logger.info("Pass 2: extracting LFCs for the selected pairs from all runs...")
        wide = pd.concat([extract_pairs(dirs[r], r, target_pairs) for r in runs], axis=1)
        wide['in_union_detected'] = wide.index.isin(union_detected)
        wide['in_random'] = wide.index.isin(rand_pairs)
        for r in runs:
            wide[f'sig_{r}'] = wide.index.isin(sig_sets[r])
        wide.to_parquet(lfc_file)
        logger.info(f"Cached LFC matrix ({wide.shape[0]:,} pairs) to {lfc_file}")

    # %%
    # Recovery of the published edges, and are the LODO sets subsets of the full run?
    sig_cols = {r: f'sig_{r}' for r in runs}
    full_edges = wide.index[wide[sig_cols[FULL]]]
    rows = []
    for c in lodo:
        c_edges = wide.index[wide[sig_cols[c]]]
        shared = len(set(c_edges) & set(full_edges))
        rows.append({
            'condition': labels[c],
            'left_out_donor': c,
            'edges_lodo': len(c_edges),
            'edges_full_run': len(full_edges),
            'shared_with_full': shared,
            'pct_of_full_recovered': round(100 * shared / len(full_edges), 1),
            'pct_of_lodo_in_full': round(100 * shared / len(c_edges), 1) if len(c_edges) else np.nan,
            'lodo_edges_not_in_full': len(c_edges) - shared,
        })
    recovery = pd.DataFrame(rows)
    recovery.to_csv(os.path.join(args.output_dir, "lodo_vs_full_recovery.tsv"), sep='\t', index=False)
    logger.info("\n" + recovery.to_string(index=False))

    # how many published edges are recovered by at least one / all LODO runs
    rec_any = wide[wide[sig_cols[FULL]]][[sig_cols[c] for c in lodo]]
    logger.info(f"Published edges recovered by >=1 LODO run: "
                f"{rec_any.any(axis=1).sum():,} ({100 * rec_any.any(axis=1).mean():.1f}%)")
    logger.info(f"Published edges recovered by all LODO runs: "
                f"{rec_any.all(axis=1).sum():,} ({100 * rec_any.all(axis=1).mean():.1f}%)")

    # %%
    # LFC agreement, per pair set
    pair_sets = {
        'full_run_significant': wide[sig_cols[FULL]].values,
        'union_detected': wide['in_union_detected'].values,
        'random': wide['in_random'].values & ~wide['in_union_detected'].values,
    }
    cor_rows = []
    for name, mask in pair_sets.items():
        sub = wide[mask]
        for a, b in itertools.combinations(runs, 2):
            x, y = sub[f'lfc_{a}'].values, sub[f'lfc_{b}'].values
            ok = np.isfinite(x) & np.isfinite(y)
            if ok.sum() < 10:
                continue
            cor_rows.append({
                'pair_set': name, 'run_a': labels[a], 'run_b': labels[b], 'n_pairs': int(ok.sum()),
                'pearson_r': round(float(stats.pearsonr(x[ok], y[ok])[0]), 3),
                'spearman_rho': round(float(stats.spearmanr(x[ok], y[ok])[0]), 3),
                'sign_concordance_pct': round(100 * float((np.sign(x[ok]) == np.sign(y[ok])).mean()), 1),
            })
    cors = pd.DataFrame(cor_rows)
    cors.to_csv(os.path.join(args.output_dir, "lodo_vs_full_correlations.tsv"), sep='\t', index=False)
    for name in pair_sets:
        logger.info(f"\nLFC agreement on '{name}':\n" +
                    cors[cors.pair_set == name].drop(columns='pair_set').to_string(index=False))

    # %%
    # Figure
    fig = plt.figure(figsize=(11, 7.8))
    outer = gridspec.GridSpec(2, 2, figure=fig, hspace=0.45, wspace=0.28, height_ratios=[1, 1])

    # panel A, recovery
    ax = fig.add_subplot(outer[0, 0])
    x = np.arange(len(lodo))
    ax.bar(x - 0.2, recovery['pct_of_full_recovered'], width=0.4, color=C_BAR,
           label='recall: % of published edges recovered')
    ax.bar(x + 0.2, recovery['pct_of_lodo_in_full'], width=0.4, color=C_ALT,
           label='precision: % of calls that are published edges')
    for xi, (a, b) in enumerate(zip(recovery['pct_of_full_recovered'], recovery['pct_of_lodo_in_full'])):
        ax.text(xi - 0.2, a, f'{a:.0f}', ha='center', va='bottom', fontsize=7)
        ax.text(xi + 0.2, b, f'{b:.0f}', ha='center', va='bottom', fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([f'without\ndonor {c}' for c in lodo])
    ax.set_ylabel('Percent')
    ax.set_ylim(0, 105)
    ax.legend(fontsize=7, frameon=False, loc='upper center')
    ax.set_title(f'A  Recall and precision vs the published 3 donor run (FDR < {FDR})',
                 fontsize=9, fontweight='bold', loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    # panel B, correlation heatmap on the published edges
    ax = fig.add_subplot(outer[0, 1])
    sub = cors[cors.pair_set == 'full_run_significant']
    m = pd.DataFrame(np.eye(len(runs)), index=[labels[r] for r in runs],
                     columns=[labels[r] for r in runs])
    for _, row in sub.iterrows():
        m.loc[row.run_a, row.run_b] = row.pearson_r
        m.loc[row.run_b, row.run_a] = row.pearson_r
    im = ax.imshow(m.values, cmap='Blues', vmin=0, vmax=1)
    ax.set_xticks(range(len(runs)))
    ax.set_xticklabels(m.columns, rotation=30, ha='right')
    ax.set_yticks(range(len(runs)))
    ax.set_yticklabels(m.index)
    for i in range(len(runs)):
        for j in range(len(runs)):
            ax.text(j, i, f'{m.values[i, j]:.2f}', ha='center', va='center', fontsize=8,
                    color='white' if m.values[i, j] > 0.55 else 'black')
    ax.set_title('B  LFC Pearson r on published edges', fontsize=9, fontweight='bold', loc='left')
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Pearson r')

    # panel C, LFC scatters against the full run, on the published edges
    inner = gridspec.GridSpecFromSubplotSpec(1, len(lodo), subplot_spec=outer[1, :], wspace=0.3)
    pub = wide[wide[sig_cols[FULL]]]
    for i, c in enumerate(lodo):
        ax = fig.add_subplot(inner[i])
        x, y = pub[f'lfc_{FULL}'].values, pub[f'lfc_{c}'].values
        ok = np.isfinite(x) & np.isfinite(y)
        ax.hexbin(x[ok], y[ok], gridsize=55, bins='log', cmap='Blues', mincnt=1)
        lim = np.nanpercentile(np.abs(np.concatenate([x[ok], y[ok]])), 99.5)
        ax.plot([-lim, lim], [-lim, lim], color='#D62728', ls='--', lw=1)
        ax.axhline(0, color=C_GREY, lw=0.5)
        ax.axvline(0, color=C_GREY, lw=0.5)
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        r = stats.pearsonr(x[ok], y[ok])[0]
        rho = stats.spearmanr(x[ok], y[ok])[0]
        conc = 100 * (np.sign(x[ok]) == np.sign(y[ok])).mean()
        ax.set_xlabel('LFC full run (3 donors)')
        ax.set_ylabel(f'LFC without donor {c}' if i == 0 else '')
        ax.set_title(f'{"C  " if i == 0 else ""}without donor {c}\n'
                     f'r={r:.2f}, rho={rho:.2f}, sign {conc:.0f}%  (n={ok.sum():,})',
                     fontsize=8, fontweight='bold' if i == 0 else 'normal', loc='left')
        for s in ['top', 'right']:
            ax.spines[s].set_visible(False)

    for ext in ['pdf', 'png']:
        out = os.path.join(args.output_dir, f"lodo_vs_full_comparison.{ext}")
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Wrote {out}")

    logger.info("All done!")


if __name__ == "__main__":
    main()
