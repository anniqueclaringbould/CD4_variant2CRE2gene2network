"""Expression-matched random-seed null for the headline "56% of genes within two hops".

The reviewer's objection is that reaching 4,356 of 7,815 genes in two hops may be a property of
any large perturbation set rather than of disease biology. This builds the matching null: seed the
same propagation with genes that are *not* CRE targets but have the same expression in
non-targeting-control cells, and ask how often a random seed set reaches as far.

Matched on what, and why:
  * expression in NTC cells (`ntc_mean_expr`, the same table Fig. 4A uses). A gene's expression
    governs how many DEGs its knockdown can produce, so an unmatched null would mostly measure
    that the CRE targets are well expressed.
  * perturbed in the promoter screen. Only 85 of the 89 CRE target genes were perturbed; the other
    four contribute no outgoing edges and are dropped from both arms, so real and null seed the
    same number of usable genes.

Everything else is held identical to b06_cascade: the same gated edge table (|log2FC| >= 0.2, FDR
0.1), the same two hops, the same `exclude_seen=False` convention, and the same denominator.

The MHC region is handled as a reported variant rather than a default. Nine of the 89 CRE target
genes lie in chr6:28-34 Mb (TNF, LTA, LTB, PSMB9, BAG6, CSNK2B, APOM, NCR3, HLA-DRB5), so it was
NOT excluded from the published figure; `--drop-mhc` re-runs both arms without them to show whether
the result rests on that one locus. MHC genes are excluded from the sampling pool either way, so a
null draw can never be dominated by it.

Outputs (to --output-dir):
    twohop_null.tsv            one row per draw
    twohop_null_summary.tsv    the real value, the null distribution, the empirical p
    twohop_null.pdf / .png
"""

from __future__ import annotations

import argparse
import os
import sys

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

DATA = "/g/stegle/schrod/code/TCell_explorer/data"
KNOCKDOWN = "/g/stegle/schrod/code/TCell/results/promoter_knockdown_efficiency.csv"

#: b06_cascade / canonical.NETWORK_LFC_MIN. Gating reproduces the published 815 / 4,356.
LFC_MIN = 0.2
HOPS = 2
#: chr6 MHC, generously bounded. Only used to define the sampling pool and the --drop-mhc arm.
MHC = ('chr6', 28_000_000, 34_000_000)
N_BINS = 10
#: `--match degseq`: how many pool genes either side of a seed's out-degree are candidates.
DEGSEQ_POOL = 25


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--output-dir', default="/g/stegle/schrod/code/TCell/Network/output")
    p.add_argument('--n-draws', type=int, default=200)
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--exclude-cre-targets', action='store_true',
                   help="drop every gene with a significant cis CRE link from the sampling pool, "
                        "not just the 89 high-confidence enhancer-like targets")
    p.add_argument('--match', choices=['expression', 'degs', 'degseq'], default='expression',
                   help="what to match the random seed sets on: NTC expression, or the number of "
                        "DEGs the seed's own knockdown produces (its out-degree in the gated "
                        "edge table)")
    p.add_argument('--enrich', choices=['none', 'Cytopus_v1.3', 'GO_Biological_Process_2025'],
                   default='none',
                   help="also score each network by how enriched its hop-1 gene set is for this "
                        "gene-set library. This asks whether random networks of the same size and "
                        "degree carry the same biology, which is a different question from how "
                        "far they reach.")
    p.add_argument('--enrich-fdr', type=float, default=0.05)
    p.add_argument('--drop-hub-top', type=int, default=0,
                   help="remove the N perturbations with the largest gated out-degree from the "
                        "network entirely -- as sources, as responses, from the seed sets, from "
                        "the sampling pool and from the denominator. Tests whether the reach is "
                        "carried by a few broadly acting regulators on both sides of the "
                        "comparison.")
    p.add_argument('--drop-mhc', action='store_true',
                   help="also drop MHC-region genes from the real seed set")
    return p.parse_args()


def load():
    """Read the built tables with pandas; every file here is a few MB."""
    deg = pd.read_parquet(f"{DATA}/deg_hits.parquet", columns=['perturbation', 'response', 'lfc'])
    # Both gene columns are dictionary-encoded with different category sets, which pandas refuses
    # to compare; plain strings are also what the adjacency dict wants.
    for col in ('perturbation', 'response'):
        deg[col] = deg[col].astype(str)
    edges = deg[(deg['lfc'].abs() >= LFC_MIN) & (deg['perturbation'] != deg['response'])]
    links = pd.read_parquet(f"{DATA}/cre_gene_links.parquet",
                            columns=['gene', 'enhancer_like_hc', 'significant_cis'])
    links['gene'] = links['gene'].astype(str)
    seeds = sorted(set(links.loc[links['enhancer_like_hc'], 'gene']))
    cre_targets = set(links.loc[links['significant_cis'], 'gene'])
    perturbed = set(pd.read_parquet(f"{DATA}/perturbation_summary.parquet",
                                    columns=['perturbation'])['perturbation'].astype(str))
    ann = pd.read_parquet(f"{DATA}/gene_annot.parquet", columns=['gene', 'chrom', 'tss'])
    mhc = set(ann.loc[(ann['chrom'].astype(str) == MHC[0])
                      & ann['tss'].between(MHC[1], MHC[2]), 'gene'].astype(str))
    expr = (pd.read_csv(KNOCKDOWN, usecols=['gene', 'ntc_mean_expr'])
            .drop_duplicates('gene').set_index('gene')['ntc_mean_expr'])
    return edges, seeds, perturbed, mhc, expr, cre_targets


def enrichment_stats(query: set[str], terms: dict[str, set[str]], universe: set[str],
                     fdr: float) -> dict:
    """How much biology a gene set carries: hypergeometric over `terms`, BH-corrected.

    Reported as the number of terms clearing `fdr` and the strength of the best one, so a network
    that is enriched for a few programmes very strongly and one that is diffusely enriched for many
    are distinguishable. The universe is the propagation universe, not all genes: a term's genes
    that were never tested cannot be reached and must not count against it.
    """
    q = query & universe
    M, N = len(universe), len(q)
    if not N:
        return {'n_sig': 0, 'best_fold': np.nan, 'best_neglog10p': 0.0, 'best_term': ''}
    names, ks, ns = [], [], []
    for t, genes in terms.items():
        k = len(q & genes)
        if k:
            names.append(t)
            ks.append(k)
            ns.append(len(genes))
    if not names:
        return {'n_sig': 0, 'best_fold': np.nan, 'best_neglog10p': 0.0, 'best_term': ''}
    ks, ns = np.array(ks), np.array(ns)
    pv = stats.hypergeom.sf(ks - 1, M, ns, N)
    order = np.argsort(pv)
    ranked = pv[order]
    # Benjamini-Hochberg
    adj = np.minimum.accumulate((ranked * len(ranked) / np.arange(1, len(ranked) + 1))[::-1])[::-1]
    best = order[0]
    return {
        'n_sig': int((adj < fdr).sum()),
        'best_fold': float((ks[best] / N) / (ns[best] / M)),
        'best_neglog10p': float(-np.log10(max(pv[best], 1e-300))),
        'best_term': names[best],
    }


def reach(adj: dict[str, set[str]], seeds: set[str], hops: int = HOPS) -> list[set[str]]:
    """Genes reached at each hop. `exclude_seen=False`, as in b06_cascade."""
    out, frontier = [], set(seeds)
    for _ in range(hops):
        nxt = set()
        for g in frontier:
            nxt |= adj.get(g, set())
        nxt.discard(None)
        out.append(set(nxt))
        frontier = nxt
        if not frontier:
            out += [set()] * (hops - len(out))
            break
    # Per-hop sets, NOT cumulative and NOT newly-reached: `exclude_seen=False` in b06_cascade, so a
    # gene reached at both hops is counted in both, and hop 2 is everything the hop-1 frontier
    # reaches. This is the definition behind the published 815 / 4,356; the union of the two hops
    # is 4,511, which is a different (and larger) number.
    return out


def main() -> int:
    a = parse_args()
    os.makedirs(a.output_dir, exist_ok=True)
    edges, seeds89, perturbed, mhc, expr, cre_targets = load()

    if a.drop_hub_top:
        deg_out = edges.groupby('perturbation')['response'].nunique().sort_values(ascending=False)
        hubs = set(deg_out.index[:a.drop_hub_top])
        print(f"dropping the {len(hubs)} largest hubs from the network "
              f"(out-degree {int(deg_out.iloc[0]):,} down to "
              f"{int(deg_out.iloc[min(a.drop_hub_top, len(deg_out)) - 1]):,}): "
              f"{', '.join(sorted(hubs)[:8])}{' ...' if len(hubs) > 8 else ''}")
        hit = sorted(hubs & set(seeds89))
        print(f"  of which CRE target genes: {len(hit)} "
              f"({', '.join(hit) if hit else 'none'})")
        edges = edges[~edges['perturbation'].isin(hubs) & ~edges['response'].isin(hubs)]
        seeds89 = [g for g in seeds89 if g not in hubs]
        perturbed = perturbed - hubs
        mhc = mhc | hubs      # keeps them out of the sampling pool too

    adj: dict[str, set[str]] = {}
    for p_, r_ in zip(edges['perturbation'].values, edges['response'].values):
        adj.setdefault(p_, set()).add(r_)

    # The denominator: every gene the propagation could possibly reach, i.e. everything that is
    # either perturbed or seen as a response in the gated edge table. Identical for both arms.
    universe = set(edges['perturbation']) | set(edges['response'])
    n_universe = len(universe)

    real_seeds = [g for g in seeds89 if g in perturbed]
    dropped_mhc = sorted(set(real_seeds) & mhc)
    if a.drop_mhc:
        real_seeds = [g for g in real_seeds if g not in mhc]
    real_seeds = [g for g in real_seeds if g in expr.index]

    print(f"universe (perturbed or responding, |log2FC| >= {LFC_MIN}): {n_universe:,}")
    print(f"CRE target genes: {len(seeds89)}  ->  perturbed: {sum(g in perturbed for g in seeds89)}"
          f"  ->  with NTC expression: {len(real_seeds)}")
    print(f"MHC-region genes among the seeds: {len(dropped_mhc)} "
          f"({', '.join(dropped_mhc) if dropped_mhc else 'none'})"
          f"{'  [DROPPED]' if a.drop_mhc else '  [kept, as in the published figure]'}")

    terms: dict[str, set[str]] = {}
    if a.enrich != 'none':
        gs = pd.read_parquet(f"{DATA}/gene_sets.parquet")
        gs = gs[gs['library'].astype(str) == a.enrich]
        for t, grp in gs.groupby(gs['term'].astype(str)):
            g = set(grp['gene'].astype(str)) & universe
            if 5 <= len(g) <= 1000:          # skip degenerate and near-universal terms
                terms[t] = g
        print(f"\nenrichment library: {a.enrich} -- {len(terms):,} terms "
              f"with 5-1000 genes inside the universe")

    real = reach(adj, set(real_seeds))
    real_pct = [100 * len(h) / n_universe for h in real]
    def score(h):
        """Enrichment of both hops. Hop 2 is 56% of the universe, so its enrichment is weak by
        construction -- it is reported so that the dilution is visible rather than assumed."""
        if not terms:
            return {}
        a1 = enrichment_stats(h[0], terms, universe, a.enrich_fdr)
        a2 = enrichment_stats(h[1], terms, universe, a.enrich_fdr)
        return {**a1, **{f'{k}_hop2': v for k, v in a2.items()}}

    real_enr = score(real) if terms else None
    print(f"\nREAL  hop1 {len(real[0]):,}  hop2 {len(real[1]):,}  "
          f"= {real_pct[1]:.1f}% of {n_universe:,}")
    if real_enr:
        print(f"      hop-1 enrichment: {real_enr['n_sig']} terms at FDR {a.enrich_fdr}; "
              f"best '{real_enr['best_term']}' "
              f"fold {real_enr['best_fold']:.2f}, -log10 p {real_enr['best_neglog10p']:.1f}")
        print(f"      hop-2 enrichment: {real_enr['n_sig_hop2']} terms at FDR {a.enrich_fdr}; "
              f"best '{real_enr['best_term_hop2']}' "
              f"fold {real_enr['best_fold_hop2']:.2f}, "
              f"-log10 p {real_enr['best_neglog10p_hop2']:.1f}")

    # ── matched sampling pool ───────────────────────────────────────────────────
    # `expression`: NTC expression, the covariate that governs how many DEGs a knockdown *can*
    #   produce. Leaves hop-1 reach free to differ, so it tests reach as a whole.
    # `degs`: the seed's own out-degree in this same gated edge table. This fixes hop 1 by
    #   construction, so any remaining difference at hop 2 is about *where in the network* the
    #   seeds sit rather than how productive their knockdowns are. The stricter test, and the one
    #   that isolates network position.
    banned = set(seeds89) | mhc | (cre_targets if a.exclude_cre_targets else set())
    # `sorted`, not set-iteration order: `perturbed` is a set, and CPython's per-process string
    # hash randomisation makes its iteration order differ between runs. Without this, --seed fixes
    # the draws but not the list they index into, so two runs of the same command sample different
    # (degree-equivalent) genes and report slightly different hop counts.
    source = sorted(expr.index) if a.match == 'expression' else sorted(perturbed)
    eligible = [g for g in source if g in perturbed and g not in banned]
    if a.match == 'degseq':
        # Nearest-degree matching: for each real seed, draw a pool gene whose out-degree is as
        # close as possible, sampling among the `DEGSEQ_POOL` nearest so the draws are not all
        # identical. Seeds are matched largest-degree first, because the hubs are the scarce
        # resource -- matching SMARCE1 (519) last would leave nothing comparable available.
        pool_deg = pd.Series({g: len(adj.get(g, ())) for g in eligible}).sort_values()
        seed_deg = sorted(((len(adj.get(g, ())), g) for g in real_seeds), reverse=True)
        print(f"\nmatching on: degree sequence")
        print(f"sampling pool: {len(pool_deg):,} genes; per draw {len(seed_deg)} genes, "
              f"each drawn from the {DEGSEQ_POOL} nearest by out-degree")
        print(f"  real seed out-degrees, top 5: "
              f"{[d for d, _ in seed_deg[:5]]}  (sum {sum(d for d, _ in seed_deg):,})")

        rng = np.random.default_rng(a.seed)
        rows = []
        order = np.array(pool_deg.values)
        names = np.array(pool_deg.index)
        for i in range(a.n_draws):
            taken, drawn = set(), []
            for d, _ in seed_deg:
                j = int(np.searchsorted(order, d))
                lo, hi = max(0, j - DEGSEQ_POOL), min(len(order), j + DEGSEQ_POOL)
                cand = [k for k in range(lo, hi) if names[k] not in taken]
                if not cand:
                    continue
                k = int(rng.choice(cand))
                taken.add(names[k])
                drawn.append(names[k])
            h = reach(adj, set(drawn))
            rows.append({'draw': i, 'n_seeds': len(drawn),
                         'n_seeds_with_degs': sum(1 for g in drawn if adj.get(g)),
                         'seed_degree_sum': int(sum(len(adj.get(g, ())) for g in drawn)),
                         'hop1': len(h[0]), 'hop2': len(h[1]),
                         'pct_hop2': 100 * len(h[1]) / n_universe,
                        **score(h)})
        null = pd.DataFrame(rows)
    elif a.match == 'expression':
        eligible = [g for g in eligible if g in expr.index]
        stat = pd.Series({g: np.log1p(expr[g]) for g in eligible})
        seed_stat = np.log1p(expr[real_seeds].values)
    else:
        stat = pd.Series({g: np.log1p(len(adj.get(g, ()))) for g in eligible})
        seed_stat = np.log1p(np.array([len(adj.get(g, ())) for g in real_seeds], dtype=float))

    if a.match != 'degseq':
        # A zero bin of its own: 44 of the 85 real seeds produce no DEG at all, and pooling them with
        # the low-but-nonzero genes would quietly give every null draw more productive seeds than the
        # real set has. With `--match degs` this reproduces "41 of 85 seeds have >=1 DEG" exactly.
        zero = stat == 0
        bin_of = pd.Series(0, index=stat.index, dtype=int)
        nz = stat[~zero]
        if len(nz):
            bin_of[nz.index] = pd.qcut(nz, min(N_BINS - 1, nz.nunique()),
                                       labels=False, duplicates='drop') + 1
        by_bin = {b: np.array(bin_of.index[bin_of == b]) for b in sorted(set(bin_of))}

        cuts = np.array([nz[bin_of[nz.index] == b].min() for b in sorted(set(bin_of[nz.index]))])
        seed_bin = np.where(seed_stat == 0, 0, np.digitize(seed_stat, cuts[1:], right=False) + 1)
        seed_bin = np.clip(seed_bin, 0, max(by_bin))
        want = pd.Series(seed_bin).value_counts().to_dict()
        print(f"\nmatching on: {a.match}")
        print(f"sampling pool: {len(stat):,} perturbed non-target non-MHC genes, "
              f"{len(by_bin)} bins ({int(zero.sum()):,} with a zero value); "
              f"per draw {sum(want.values())} genes")

        rng = np.random.default_rng(a.seed)
        rows = []
        for i in range(a.n_draws):
            drawn: list[str] = []
            for b, n in want.items():
                cand = by_bin[b]
                drawn += list(rng.choice(cand, size=min(n, len(cand)), replace=False))
            h = reach(adj, set(drawn))
            rows.append({'draw': i, 'n_seeds': len(drawn),
                         'n_seeds_with_degs': sum(1 for g in drawn if adj.get(g)),
                         'seed_degree_sum': int(sum(len(adj.get(g, ())) for g in drawn)),
                         'hop1': len(h[0]), 'hop2': len(h[1]),
                         'pct_hop2': 100 * len(h[1]) / n_universe,
                        **score(h)})
        null = pd.DataFrame(rows)
    null.to_csv(os.path.join(a.output_dir, "twohop_null.tsv"), sep='\t', index=False)

    ge = int((null['pct_hop2'] >= real_pct[1]).sum())
    p_emp = (ge + 1) / (len(null) + 1)
    print(f"\nNULL  seeds with >=1 DEG: median {null['n_seeds_with_degs'].median():.0f}"
          f" of {len(real_seeds)}  (real: {sum(1 for g in real_seeds if adj.get(g))})")
    print(f"      seed out-degree sum: median {null['seed_degree_sum'].median():,.0f}"
          f"  (real: {sum(len(adj.get(g, ())) for g in real_seeds):,})")
    print(f"      hop1 reach: median {null['hop1'].median():,.0f}  (real: {len(real[0]):,})")
    print(f"NULL  hop2 reach: median {null['pct_hop2'].median():.1f}%  "
          f"mean {null['pct_hop2'].mean():.1f}%  "
          f"range {null['pct_hop2'].min():.1f}-{null['pct_hop2'].max():.1f}%")
    if real_enr:
        ge_e = int((null['n_sig'] >= real_enr['n_sig']).sum())
        print(f"\nENRICHMENT  hop-1 terms at FDR {a.enrich_fdr}: real {real_enr['n_sig']}  |  "
              f"null median {null['n_sig'].median():.0f}, "
              f"range {null['n_sig'].min():.0f}-{null['n_sig'].max():.0f}")
        print(f"            draws with >= the real count: {ge_e}/{len(null)}  "
              f"(empirical p = {(ge_e + 1) / (len(null) + 1):.3f})")
        ge2 = int((null['n_sig_hop2'] >= real_enr['n_sig_hop2']).sum())
        print(f"            hop-2 terms: real {real_enr['n_sig_hop2']}  |  "
              f"null median {null['n_sig_hop2'].median():.0f}, "
              f"range {null['n_sig_hop2'].min():.0f}-{null['n_sig_hop2'].max():.0f}  |  "
              f"{ge2}/{len(null)} draws >= real "
              f"(empirical p = {(ge2 + 1) / (len(null) + 1):.3f})")
    print(f"      draws reaching >= the real {real_pct[1]:.1f}%: {ge}/{len(null)}  "
          f"(empirical p = {p_emp:.3f})")

    pd.DataFrame([{
        'n_universe': n_universe, 'n_seeds': len(real_seeds),
        'match_on': a.match,
        'mhc_in_seeds': len(dropped_mhc), 'mhc_dropped': a.drop_mhc,
        'real_hop1': len(real[0]), 'real_hop2': len(real[1]), 'real_pct_hop2': real_pct[1],
        'null_median_pct': null['pct_hop2'].median(), 'null_mean_pct': null['pct_hop2'].mean(),
        'null_sd_pct': null['pct_hop2'].std(), 'n_draws': len(null),
        'n_ge_real': ge, 'p_empirical': p_emp,
    }]).to_csv(os.path.join(a.output_dir, "twohop_null_summary.tsv"), sep='\t', index=False)

    # ── figure ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(4.4, 2.6))
    ax.hist(null['pct_hop2'], bins=25, color='#B0B0B0', edgecolor='white', linewidth=0.4)
    ax.axvline(real_pct[1], color='#D62728', lw=1.8)
    ax.annotate(f"CRE targets\n{real_pct[1]:.0f}%", (real_pct[1], ax.get_ylim()[1]),
                textcoords='offset points', xytext=(5, -8), ha='left', va='top',
                fontsize=7.5, color='#D62728', fontweight='bold')
    ax.set_xlabel(f'genes reached within {HOPS} hops (% of {n_universe:,})')
    ax.set_ylabel(f'random seed sets (n={len(null)})')
    ax.set_title({'expression': 'matched on NTC expression', 'degs': 'matched on DEG count',
                  'degseq': 'matched on degree sequence'}[a.match],
                 fontsize=8, loc='left')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    fig.tight_layout()
    for ext in ['pdf', 'png']:
        out = os.path.join(a.output_dir, f"twohop_null.{ext}")
        fig.savefig(out, dpi=300, bbox_inches='tight')
        print(f"Wrote {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
