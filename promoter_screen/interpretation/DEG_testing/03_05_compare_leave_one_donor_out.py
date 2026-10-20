# %%
"""
Compare the differentially expressed edges detected in the leave-one-donor-out (LODO) runs of
03_01_run_DESeq2_batched.py.

An edge is a (perturbation, response gene) pair called significant at adj_p_value < 0.1.

The number of targets that pass the --min-bulks filter differs between the LODO conditions, so the
comparison is restricted to the perturbations and genes that were tested in *all* conditions. A
perturbation that was dropped in one condition for lack of pseudobulks is excluded everywhere, since
otherwise differences in testing coverage would show up as apparent disagreement.

Outputs (to --output-dir):
    lodo_edges_near_significant.parquet   cached edges with adj_p < CACHE_P (intermediate)
    lodo_tested_contrasts.csv             which perturbation was tested in which condition
    lodo_summary.tsv                      per condition edge counts
    lodo_jaccard.tsv                      pairwise Jaccard similarity
    lodo_intersections.tsv                UpSet intersection sizes
    lodo_comparison.pdf / .png            UpSet, Jaccard, replication and near-miss panels
"""

import os
import sys
import argparse
import itertools
from functools import reduce

import numpy as np
import pandas as pd
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

# adj_p threshold defining a detected edge
FDR = 0.1
# rows below this adj_p are cached, so that near misses can be inspected without re-reading 20 GB
CACHE_P = 0.5

C_BAR = '#2C7FB8'
C_DOT = '#253494'
C_GREY = '#B0B0B0'


def load_condition(deg_dir, condition, cache_p=CACHE_P):
    """Stream the chunk files of one condition, return near significant edges and tested contrasts."""
    files = sorted(f for f in os.listdir(deg_dir) if f.endswith('.csv'))
    logger.info(f"  {condition}: reading {len(files)} chunk files from {deg_dir}")

    edges, contrasts, genes = [], set(), set()
    n_tests = 0

    for i, file_name in enumerate(files, start=1):
        df = pd.read_csv(os.path.join(deg_dir, file_name),
                         usecols=['contrast', 'variable', 'log_fc', 'adj_p_value'])
        n_tests += len(df)
        contrasts.update(df['contrast'].unique())
        genes.update(df['variable'].unique())
        edges.append(df[df['adj_p_value'] < cache_p])
        if i % 25 == 0 or i == len(files):
            logger.info(f"    {i}/{len(files)} files, {n_tests:,} tests so far")

    edges = pd.concat(edges, ignore_index=True)
    edges['condition'] = condition
    logger.info(f"  {condition}: {n_tests:,} tests, {len(contrasts):,} contrasts, "
                f"{len(genes):,} genes, {len(edges):,} rows with adj_p < {cache_p}")
    return edges, contrasts, genes, n_tests


def jaccard(a, b):
    if not a and not b:
        return np.nan
    return len(a & b) / len(a | b)


def plot_upset(ax_bar, ax_matrix, ax_sets, edge_sets, conditions):
    """Minimal UpSet plot: intersection bars, membership matrix and set size bars."""
    # all non empty exclusive intersections
    combos = []
    for r in range(len(conditions), 0, -1):
        for combo in itertools.combinations(conditions, r):
            inside = [edge_sets[c] for c in combo]
            outside = [edge_sets[c] for c in conditions if c not in combo]
            exclusive = reduce(set.intersection, inside)
            for o in outside:
                exclusive = exclusive - o
            if exclusive:
                combos.append((combo, len(exclusive)))
    combos.sort(key=lambda x: -x[1])

    x = np.arange(len(combos))
    sizes = [n for _, n in combos]

    # intersection size bars
    ax_bar.bar(x, sizes, color=C_BAR, width=0.6)
    for xi, n in zip(x, sizes):
        ax_bar.text(xi, n, f'{n:,}', ha='center', va='bottom', fontsize=7)
    ax_bar.set_ylabel('Edges in intersection')
    ax_bar.set_xlim(-0.6, len(combos) - 0.4)
    ax_bar.set_ylim(0, max(sizes) * 1.18)
    ax_bar.set_xticks([])
    for side in ['top', 'right']:
        ax_bar.spines[side].set_visible(False)

    # membership matrix
    for yi, cond in enumerate(conditions):
        ax_matrix.scatter(x, [yi] * len(x), s=45, color=C_GREY, zorder=1)
    for xi, (combo, _) in enumerate(combos):
        idx = [conditions.index(c) for c in combo]
        ax_matrix.plot([xi, xi], [min(idx), max(idx)], color=C_DOT, lw=1.6, zorder=2)
        ax_matrix.scatter([xi] * len(idx), idx, s=45, color=C_DOT, zorder=3)
    ax_matrix.set_xlim(-0.6, len(combos) - 0.4)
    ax_matrix.set_ylim(-0.6, len(conditions) - 0.4)
    ax_matrix.set_yticks(range(len(conditions)))
    ax_matrix.set_yticklabels(conditions)
    ax_matrix.set_xticks([])
    ax_matrix.invert_yaxis()
    for side in ['top', 'right', 'bottom', 'left']:
        ax_matrix.spines[side].set_visible(False)

    # set size bars, aligned with the matrix rows
    totals = [len(edge_sets[c]) for c in conditions]
    y = np.arange(len(conditions))
    ax_sets.barh(y, totals, color=C_BAR, height=0.6)
    for yi, n in zip(y, totals):
        ax_sets.text(n, yi, f'{n:,} ', ha='right', va='center', fontsize=7, color='white')
    ax_sets.set_ylim(-0.6, len(conditions) - 0.4)
    ax_sets.invert_yaxis()
    ax_sets.invert_xaxis()
    ax_sets.set_yticks([])
    ax_sets.set_xlabel('Total edges')
    for side in ['top', 'right', 'left']:
        ax_sets.spines[side].set_visible(False)

    return combos


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--input-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor",
                   help='Directory containing the pseudobulk_deseq2_without_donor_* result directories')
    p.add_argument('--output-dir', type=str, default="/g/stegle/schrod/code/TCell/Donor/output",
                   help='Directory for the summary tables and figures')
    p.add_argument('--keep-self-edges', action='store_true',
                   help='Keep on-target (perturbation == response gene) edges. They are excluded by '
                        'default, as they measure the knockdown itself rather than a regulatory edge')
    p.add_argument('--force-reload', action='store_true',
                   help='Ignore the cached edge table and re-read all chunk files')
    return p.parse_args()


def main():
    args = parse_args()
    input_dir, output_dir = args.input_dir, args.output_dir
    keep_self_edges, force_reload = args.keep_self_edges, args.force_reload

    os.makedirs(output_dir, exist_ok=True)
    cache_file = os.path.join(output_dir, "lodo_edges_near_significant.parquet")
    tested_file = os.path.join(output_dir, "lodo_tested_contrasts.csv")
    genes_file = os.path.join(output_dir, "lodo_tested_genes.csv")
    meta_file = os.path.join(output_dir, "lodo_condition_meta.csv")

    # %%
    # Find the LODO result directories
    prefix = "pseudobulk_deseq2_without_donor_"
    conditions = sorted(d.replace(prefix, "") for d in os.listdir(input_dir) if d.startswith(prefix))
    if len(conditions) < 2:
        logger.error(f"Found {len(conditions)} LODO directories in {input_dir}, need at least 2.")
        sys.exit(1)
    labels = {c: f"without donor {c}" for c in conditions}
    logger.info(f"Found LODO conditions: {', '.join(conditions)}")

    # %%
    # Load, using the cached edge table when available
    cached = all(os.path.exists(f) for f in [cache_file, tested_file, genes_file, meta_file])
    if cached and not force_reload:
        logger.info(f"Loading cached edges from {cache_file}")
        edges = pd.read_parquet(cache_file)
        tested = pd.read_csv(tested_file)
        tested_genes = pd.read_csv(genes_file)
        meta = pd.read_csv(meta_file).set_index('condition')
        tested_contrasts = {c: set(tested.loc[tested[c], 'contrast']) for c in conditions}
        gene_sets = {c: set(tested_genes.loc[tested_genes[c], 'gene']) for c in conditions}
        n_tests = {c: int(meta.loc[c, 'n_tests']) for c in conditions}
    else:
        logger.info(f"Reading DESeq2 results (this streams ~20 GB, keeping only adj_p < {CACHE_P})...")
        all_edges, tested_contrasts, gene_sets, n_tests = [], {}, {}, {}
        for c in conditions:
            e, contrasts, genes, nt = load_condition(os.path.join(input_dir, prefix + c), c)
            all_edges.append(e)
            tested_contrasts[c] = contrasts
            gene_sets[c] = genes
            n_tests[c] = nt
        edges = pd.concat(all_edges, ignore_index=True)

        # cache which contrast and which gene was tested in which condition
        tested = pd.DataFrame({'contrast': sorted(set().union(*tested_contrasts.values()))})
        for c in conditions:
            tested[c] = tested['contrast'].isin(tested_contrasts[c])
        tested.to_csv(tested_file, index=False)

        tested_genes = pd.DataFrame({'gene': sorted(set().union(*gene_sets.values()))})
        for c in conditions:
            tested_genes[c] = tested_genes['gene'].isin(gene_sets[c])
        tested_genes.to_csv(genes_file, index=False)

        pd.DataFrame({'condition': conditions,
                      'n_tests': [n_tests[c] for c in conditions],
                      'n_contrasts': [len(tested_contrasts[c]) for c in conditions],
                      'n_genes': [len(gene_sets[c]) for c in conditions]}).to_csv(meta_file, index=False)

        edges.to_parquet(cache_file, index=False)
        logger.info(f"Cached {len(edges):,} near significant edges to {cache_file}")

    # %%
    # Restrict to what was tested everywhere. A perturbation missing in one condition (too few
    # pseudobulks) is dropped from all conditions.
    common_contrasts = reduce(set.intersection, tested_contrasts.values())
    common_genes = reduce(set.intersection, gene_sets.values())
    logger.info("Tested perturbations per condition: " +
                ", ".join(f"{labels[c]}={len(tested_contrasts[c]):,}" for c in conditions))
    logger.info(f"Perturbations tested in all conditions: {len(common_contrasts):,}")
    for c in conditions:
        dropped = len(tested_contrasts[c] - common_contrasts)
        logger.info(f"  {labels[c]}: {dropped:,} perturbations tested here but not everywhere")

    edges = edges[edges['contrast'].isin(common_contrasts) & edges['variable'].isin(common_genes)]

    if not keep_self_edges:
        n_self = (edges['contrast'] == edges['variable']).sum()
        edges = edges[edges['contrast'] != edges['variable']]
        logger.info(f"Excluded {n_self:,} on-target rows (perturbation == response gene)")

    # %%
    # Detected edges per condition
    sig = edges[edges['adj_p_value'] < FDR].copy()
    sig['edge'] = sig['contrast'] + '->' + sig['variable']
    edge_sets = {c: set(sig.loc[sig.condition == c, 'edge']) for c in conditions}

    union = reduce(set.union, edge_sets.values())
    core = reduce(set.intersection, edge_sets.values())

    summary = pd.DataFrame({
        'condition': [labels[c] for c in conditions],
        'left_out_donor': conditions,
        'perturbations_tested': [len(tested_contrasts[c]) for c in conditions],
        'perturbations_tested_in_all': [len(common_contrasts)] * len(conditions),
        'tests_performed': [n_tests[c] for c in conditions],
        'edges_detected': [len(edge_sets[c]) for c in conditions],
        'edges_unique_to_condition': [
            len(edge_sets[c] - reduce(set.union, [edge_sets[o] for o in conditions if o != c]))
            for c in conditions],
        'perturbations_with_edge': [sig.loc[sig.condition == c, 'contrast'].nunique() for c in conditions],
    })
    summary['pct_of_union'] = (100 * summary['edges_detected'] / len(union)).round(1)
    summary.to_csv(os.path.join(output_dir, "lodo_summary.tsv"), sep='\t', index=False)
    logger.info("\n" + summary.to_string(index=False))
    logger.info(f"Union of detected edges: {len(union):,}, detected in all conditions: {len(core):,} "
                f"({100 * len(core) / len(union):.1f}% of union)")

    # %%
    # Pairwise Jaccard similarity
    jac = pd.DataFrame(index=conditions, columns=conditions, dtype=float)
    for a, b in itertools.product(conditions, conditions):
        jac.loc[a, b] = 1.0 if a == b else jaccard(edge_sets[a], edge_sets[b])
    jac.to_csv(os.path.join(output_dir, "lodo_jaccard.tsv"), sep='\t')
    offdiag = [jac.loc[a, b] for a, b in itertools.combinations(conditions, 2)]
    logger.info(f"Pairwise Jaccard: " +
                ", ".join(f"{a} vs {b} = {jac.loc[a, b]:.3f}" for a, b in itertools.combinations(conditions, 2)))
    logger.info(f"Mean pairwise Jaccard: {np.mean(offdiag):.3f}")

    # %%
    # Is the low overlap driven by weak, near threshold perturbations? Stratify by how many edges a
    # perturbation has, keeping only perturbations that reach the threshold in every condition.
    n_edges_per_pert = sig.groupby(['condition', 'contrast']).size().unstack(0) \
                          .reindex(sorted(common_contrasts)).fillna(0)
    rows = []
    for thr in [0, 5, 10, 25, 50]:
        strong = set(n_edges_per_pert.index[(n_edges_per_pert >= thr).all(axis=1)])
        sub = {c: {e for e in edge_sets[c] if e.split('->')[0] in strong} for c in conditions}
        js = [jaccard(sub[a], sub[b]) for a, b in itertools.combinations(conditions, 2)]
        u, k = reduce(set.union, sub.values()), reduce(set.intersection, sub.values())
        rows.append({'min_edges_per_perturbation_in_all_conditions': thr,
                     'n_perturbations': len(strong), 'mean_jaccard': round(float(np.mean(js)), 3),
                     'edges_union': len(u), 'edges_in_all_conditions': len(k),
                     'pct_in_all_conditions': round(100 * len(k) / len(u), 1) if u else np.nan})
    strength = pd.DataFrame(rows)
    strength.to_csv(os.path.join(output_dir, "lodo_jaccard_by_strength.tsv"), sep='\t', index=False)
    logger.info("Jaccard stratified by perturbation strength:\n" + strength.to_string(index=False))

    # %%
    # How often is an edge replicated, and are the unique ones just below threshold?
    n_cond = pd.Series({e: sum(e in edge_sets[c] for c in conditions) for e in union})
    replication = n_cond.value_counts().sort_index()

    # sign concordance of the log fold change among edges detected in all conditions
    lfc_wide = sig.pivot_table(index='edge', columns='condition', values='log_fc')
    lfc_core = lfc_wide.loc[sorted(core)] if core else lfc_wide.iloc[:0]
    if len(lfc_core):
        signs = np.sign(lfc_core[conditions].values)
        concordant = (np.abs(signs.sum(axis=1)) == len(conditions)).mean()
        logger.info(f"Sign concordance among edges detected in all conditions: {100 * concordant:.2f}%")
    else:
        concordant = np.nan

    # for edges detected in exactly one condition, adj_p in the conditions that missed them
    adjp_wide = edges.assign(edge=edges['contrast'] + '->' + edges['variable']) \
                     .pivot_table(index='edge', columns='condition', values='adj_p_value')
    unique_edges = set(n_cond[n_cond == 1].index)
    missed_p = []
    for c in conditions:
        uniq_c = sorted(edge_sets[c] & unique_edges)
        others = [o for o in conditions if o != c]
        vals = adjp_wide.reindex(uniq_c)[others].values.ravel()
        # rows absent from the cache had adj_p >= CACHE_P, treat them as such
        vals = np.where(np.isnan(vals), CACHE_P, vals)
        missed_p.append(vals)
    missed_p = np.concatenate(missed_p) if missed_p else np.array([])

    inter_tbl = pd.DataFrame({'n_conditions': replication.index, 'n_edges': replication.values})
    inter_tbl['pct_of_union'] = (100 * inter_tbl['n_edges'] / len(union)).round(1)
    inter_tbl.to_csv(os.path.join(output_dir, "lodo_replication.tsv"), sep='\t', index=False)
    logger.info("\n" + inter_tbl.to_string(index=False))

    # %%
    # Figure
    fig = plt.figure(figsize=(11, 7.5))
    outer = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.32,
                              height_ratios=[1.15, 1])

    # panel A, UpSet
    inner = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=outer[0, :], hspace=0.06, wspace=0.04,
                                             height_ratios=[3, 1], width_ratios=[1, 4])
    ax_bar = fig.add_subplot(inner[0, 1])
    ax_matrix = fig.add_subplot(inner[1, 1], sharex=ax_bar)
    ax_sets = fig.add_subplot(inner[1, 0], sharey=ax_matrix)
    combos = plot_upset(ax_bar, ax_matrix, ax_sets, edge_sets, conditions)
    ax_bar.set_title(f'A  Differentially expressed edges per leave-one-donor-out condition '
                     f'(FDR < {FDR}, {len(common_contrasts):,} commonly tested perturbations)',
                     fontsize=9, fontweight='bold', loc='left')

    inter_out = pd.DataFrame({
        'intersection': [' & '.join(labels[c] for c in combo) for combo, _ in combos],
        'n_edges': [n for _, n in combos]})
    inter_out.to_csv(os.path.join(output_dir, "lodo_intersections.tsv"), sep='\t', index=False)

    # panel B, Jaccard heatmap
    ax = fig.add_subplot(outer[1, 0])
    m = jac.loc[conditions, conditions].values.astype(float)
    im = ax.imshow(m, cmap='Blues', vmin=0, vmax=1)
    ax.set_xticks(range(len(conditions)))
    ax.set_xticklabels([labels[c] for c in conditions], rotation=30, ha='right')
    ax.set_yticks(range(len(conditions)))
    ax.set_yticklabels([labels[c] for c in conditions])
    for i in range(len(conditions)):
        for j in range(len(conditions)):
            ax.text(j, i, f'{m[i, j]:.2f}', ha='center', va='center', fontsize=9,
                    color='white' if m[i, j] > 0.55 else 'black')
    ax.set_title(f'B  Jaccard similarity of edge sets\n(mean off-diagonal {np.mean(offdiag):.2f})',
                 fontsize=9, fontweight='bold', loc='left')
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label='Jaccard index')

    # panel C, replication and near miss
    inner2 = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1, 1], wspace=0.45)

    ax = fig.add_subplot(inner2[0])
    ax.bar(inter_tbl['n_conditions'], inter_tbl['n_edges'], color=C_BAR, width=0.6)
    for xi, n, p in zip(inter_tbl['n_conditions'], inter_tbl['n_edges'], inter_tbl['pct_of_union']):
        ax.text(xi, n, f'{p:.0f}%', ha='center', va='bottom', fontsize=7)
    ax.set_xticks(inter_tbl['n_conditions'])
    ax.set_xlabel('Conditions detecting the edge')
    ax.set_ylabel('Edges')
    ax.set_ylim(0, inter_tbl['n_edges'].max() * 1.15)
    ax.set_title('C  Replication', fontsize=9, fontweight='bold', loc='left')
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)

    ax = fig.add_subplot(inner2[1])
    if len(missed_p):
        ax.hist(missed_p, bins=np.linspace(0, CACHE_P, 26), color=C_GREY, edgecolor='white')
        ax.axvline(FDR, color='#D62728', ls='--', lw=1.2)
        ax.text(FDR, ax.get_ylim()[1] * 0.95, f' FDR {FDR}', color='#D62728', fontsize=7,
                va='top', ha='left')
        below = (missed_p < 0.25).mean()
        ax.set_title(f'D  adj. p of condition-unique edges\nin the conditions that missed them\n'
                     f'({100 * below:.0f}% still below 0.25)',
                     fontsize=9, fontweight='bold', loc='left')
    ax.set_xlabel(f'adj. p value (capped at {CACHE_P})')
    ax.set_ylabel('Edges')
    for side in ['top', 'right']:
        ax.spines[side].set_visible(False)

    for ext in ['pdf', 'png']:
        out = os.path.join(output_dir, f"lodo_comparison.{ext}")
        fig.savefig(out, dpi=300, bbox_inches='tight')
        logger.info(f"Wrote {out}")

    logger.info("All done!")


if __name__ == "__main__":
    main()
