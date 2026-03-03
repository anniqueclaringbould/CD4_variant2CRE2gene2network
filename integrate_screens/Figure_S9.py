#%%
"""
Panel A standalone — CRE functional annotation rate (Cytopus).

Reads TAP-seq CRE screen results directly from Table S5.xlsx and
computational predictions from Table S2.xlsx.

Three cascade depths are compared (L1 only, +hop 1, +hop 2),
each with real Perturb-seq, random-gene null, and shuffled-cascade null.

Input files (in ./input/):
  - Table S5.xlsx              — TAP-seq CRE perturbation screen results
  - Table S2.xlsx              — computational CRE-gene predictions
  - gene_set.csv               — variable gene set (background)
  - perturbseq_hits_nolfc.csv  — Perturb-seq gene-gene hits
  - geneset_cache/             — cached Cytopus gene sets

Output (in ./output/):
  - panel_A_standalone.{png,pdf,svg}
"""

import pandas as pd
import numpy as np
import gseapy as gp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from pathlib import Path
import warnings, json, logging

warnings.filterwarnings("ignore")
logging.getLogger("gseapy").setLevel(logging.CRITICAL)

plt.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
})

INPUT_DIR = OUTPUT_DIR = Path("/g/stegle/schrod/code/TCell/Figure_S9")

# ── Parameters ──
TOP_N     = 20
TOP_N_VIZ = 10
LFC_CUTOFF = 0.2
SIG_CUTOFF = 0.05
N_PERM     = 100


# ═══════════════════════════════════════════════════════════════════
#  Data loading
# ═══════════════════════════════════════════════════════════════════

def load_geneset_library_cached(lib_name):
    cache_dir = INPUT_DIR / "geneset_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_file = cache_dir / f"{lib_name}.json"
    if cache_file.exists():
        with open(cache_file) as f:
            return json.load(f)
    lib = gp.get_library(lib_name)
    with open(cache_file, "w") as f:
        json.dump({k: list(v) if isinstance(v, set) else v
                   for k, v in lib.items()}, f)
    return lib


def load_data():
    print("Loading data ...")

    # Variable genes
    var_gene_set = set(
        pd.read_csv(INPUT_DIR / "gene_set.csv")["bpcells_name"].unique()
    )

    # TAP-seq from Table S5
    tapseq = pd.read_excel(
        INPUT_DIR / "Table S5.xlsx",
        sheet_name="cre_perturbation_screen_results",
    )
    enh = tapseq[
        (tapseq["enhancer_like_interaction"] == True)
        & (tapseq["significant"] == True)
    ].copy()
    enh_target_genes = set(enh["response_id"].unique())

    # Perturb-seq hits
    ps = pd.read_excel("/g/stegle/schrod/code/TCell/Table S12.xlsx")
    ps.rename(columns={"KO": "contrast", "gene_name": "variable", "adj_pval": "adj_p_value", "LFC": "log_fc"}, inplace=True)
    ps = ps[ps["contrast"] != ps["variable"]]
    ps = ps[ps["log_fc"].abs() >= LFC_CUTOFF]
    ps = ps[ps["contrast"].isin(var_gene_set) & ps["variable"].isin(var_gene_set)]
    ps = ps.sort_values("adj_p_value")
    ps["rank"] = ps.groupby("contrast").cumcount() + 1
    ps = ps[ps["rank"] <= TOP_N]

    ps_lookup = {}
    for _, row in ps.iterrows():
        ps_lookup.setdefault(row["contrast"], []).append({
            "target": row["variable"],
            "adj_p":  float(row["adj_p_value"]),
            "log_fc": float(row["log_fc"]),
        })

    # S2 computational predictions
    s2 = pd.read_excel(INPUT_DIR / "Table S2.xlsx")
    s2_lookup = {}
    for _, row in s2.iterrows():
        if row["gene_name"] in var_gene_set:
            s2_lookup.setdefault(row["peak"], set()).add(row["gene_name"])

    print(f"  {len(var_gene_set)} variable genes, "
          f"{len(enh['grna_target'].unique())} CREs, "
          f"{len(enh_target_genes)} target genes, "
          f"{len(s2_lookup)} S2 peaks")

    return enh, ps_lookup, var_gene_set, enh_target_genes, s2_lookup


# ═══════════════════════════════════════════════════════════════════
#  Cascade builders
# ═══════════════════════════════════════════════════════════════════

def _expand_hop(seeds, ps_lookup, var_gene_set, gate=None, cap=None):
    out = set()
    for gene in seeds:
        if gate is not None and gene not in gate:
            continue
        if gene not in ps_lookup:
            continue
        count = 0
        for hit in ps_lookup[gene]:
            if cap and count >= cap:
                break
            t = hit["target"]
            if var_gene_set is None or t in var_gene_set:
                out.add(t)
                count += 1
    return out


def build_tapseq_layers(cre, enh, ps_lookup, var_gene_set, enh_target_genes):
    l1 = set(enh.loc[enh["grna_target"] == cre, "response_id"])
    l2 = _expand_hop(l1, ps_lookup, None, gate=enh_target_genes)
    l3 = _expand_hop(l2, ps_lookup, var_gene_set, cap=TOP_N_VIZ)
    return l1, l2, l3


def build_s2_layers(cre, s2_lookup, ps_lookup, var_gene_set):
    l1 = s2_lookup.get(cre, set()).copy()
    l2 = _expand_hop(l1, ps_lookup, None)
    l3 = _expand_hop(l2, ps_lookup, var_gene_set, cap=TOP_N_VIZ)
    return l1, l2, l3


def build_shuffled_cascade(l1, ps_lookup, var_gene_set,
                           enh_target_genes, all_ps_sources, rng):
    sources = list(all_ps_sources)
    if enh_target_genes is not None:
        n = sum(1 for g in l1 if g in enh_target_genes and g in ps_lookup)
    else:
        n = sum(1 for g in l1 if g in ps_lookup)

    fake_l2 = set()
    for _ in range(n):
        donor = sources[rng.integers(len(sources))]
        if donor in ps_lookup:
            for hit in ps_lookup[donor]:
                fake_l2.add(hit["target"])

    fake_l3 = set()
    for gene in fake_l2:
        if gene not in ps_lookup:
            continue
        donor = sources[rng.integers(len(sources))]
        if donor in ps_lookup:
            count = 0
            for hit in ps_lookup[donor]:
                if count >= TOP_N_VIZ:
                    break
                if hit["target"] in var_gene_set:
                    fake_l3.add(hit["target"])
                    count += 1
    return fake_l2, fake_l3


# ═══════════════════════════════════════════════════════════════════
#  Enrichment helpers
# ═══════════════════════════════════════════════════════════════════

def run_hypergeom(gene_set, lib, bg):
    gl = sorted(gene_set)
    if len(gl) < 2:
        return {}
    try:
        enr = gp.enrich(gene_list=gl, gene_sets=lib, background=bg,
                        outdir=None, no_plot=True, verbose=False)
        return {r["Term"]: r["Adjusted P-value"]
                for _, r in enr.results.iterrows()}
    except Exception:
        return {}


# ═══════════════════════════════════════════════════════════════════
#  Main enrichment loop
# ═══════════════════════════════════════════════════════════════════

def run_enrichment(enh, ps_lookup, var_gene_set, enh_target_genes,
                   s2_lookup, lib, bg):
    all_cres = sorted(enh["grna_target"].unique())
    n_cres = len(all_cres)
    all_ps_sources = sorted(ps_lookup.keys())
    rng = np.random.default_rng(42)

    print(f"  {n_cres} CREs, {len(all_ps_sources)} PS sources, {N_PERM} perms\n")

    # Storage: 6 conditions x 3 depths, tracking frac_annotated
    D = range(3)
    real_nsig = {d: [] for d in D}
    real_sizes = {d: [] for d in D}
    null1_frac = {d: [] for d in D}
    null2_frac = {d: [] for d in D}
    s2_nsig = {d: [] for d in D}
    s2_frac = {d: [] for d in D}
    s2_sizes = {d: [] for d in D}
    s2_null1_frac = {d: [] for d in D}
    s2_null2_frac = {d: [] for d in D}

    # Track unique significant terms (union across all CREs)
    tap_ps_terms   = {d: set() for d in D}
    tap_rnd_terms  = {d: set() for d in D}
    tap_shuf_terms = {d: set() for d in D}
    s2_ps_terms    = {d: set() for d in D}

    for idx, cre in enumerate(all_cres):
        # ── TAP-seq cascade ──
        l1, l2, l3 = build_tapseq_layers(
            cre, enh, ps_lookup, var_gene_set, enh_target_genes)
        sets = [l1, l1 | l2, l1 | l2 | l3]
        sizes = [len(s) for s in sets]
        for d in D:
            real_sizes[d].append(sizes[d])

        # Real enrichment
        for d in D:
            res = run_hypergeom(sets[d], lib, bg)
            ns = sum(1 for p in res.values() if p < SIG_CUTOFF)
            real_nsig[d].append(ns)
            tap_ps_terms[d].update(t for t, p in res.items() if p < SIG_CUTOFF)

        # Null 1: random genes
        for d in D:
            if d == 0 or sizes[d] <= sizes[0]:
                null1_frac[d].append(1.0 if real_nsig[0][idx] > 0 else 0.0)
                continue
            n_add = sizes[d] - sizes[0]
            pool = sorted(var_gene_set - l1)
            hits = 0
            for _ in range(N_PERM):
                ri = rng.choice(len(pool), min(n_add, len(pool)), replace=False)
                pset = l1 | {pool[i] for i in ri}
                res = run_hypergeom(pset, lib, bg)
                if any(p < SIG_CUTOFF for p in res.values()):
                    hits += 1
                tap_rnd_terms[d].update(t for t, p in res.items() if p < SIG_CUTOFF)
            null1_frac[d].append(hits / N_PERM)

        # Null 2: shuffled cascade
        for d in D:
            if d == 0 or sizes[d] <= sizes[0]:
                null2_frac[d].append(1.0 if real_nsig[0][idx] > 0 else 0.0)
                continue
            hits = 0
            for _ in range(N_PERM):
                fl2, fl3 = build_shuffled_cascade(
                    l1, ps_lookup, var_gene_set, enh_target_genes,
                    all_ps_sources, rng)
                pset = (l1 | fl2) if d == 1 else (l1 | fl2 | fl3)
                res = run_hypergeom(pset, lib, bg)
                if any(p < SIG_CUTOFF for p in res.values()):
                    hits += 1
                tap_shuf_terms[d].update(t for t, p in res.items() if p < SIG_CUTOFF)
            null2_frac[d].append(hits / N_PERM)

        # ── S2 cascade ──
        sl1, sl2, sl3 = build_s2_layers(cre, s2_lookup, ps_lookup, var_gene_set)
        s2sets = [sl1, sl1 | sl2, sl1 | sl2 | sl3]
        s2sz = [len(s) for s in s2sets]
        for d in D:
            s2_sizes[d].append(s2sz[d])

        for d in D:
            res = run_hypergeom(s2sets[d], lib, bg)
            ns = sum(1 for p in res.values() if p < SIG_CUTOFF)
            s2_nsig[d].append(ns)
            s2_frac[d].append(1.0 if ns > 0 else 0.0)
            s2_ps_terms[d].update(t for t, p in res.items() if p < SIG_CUTOFF)

        # S2 null 1: random
        for d in D:
            if d == 0 or s2sz[d] <= s2sz[0]:
                s2_null1_frac[d].append(s2_frac[d][-1])
                continue
            n_add = s2sz[d] - s2sz[0]
            pool = sorted(var_gene_set - sl1)
            hits = 0
            for _ in range(N_PERM):
                ri = rng.choice(len(pool), min(n_add, len(pool)), replace=False)
                pset = sl1 | {pool[i] for i in ri}
                res = run_hypergeom(pset, lib, bg)
                if any(p < SIG_CUTOFF for p in res.values()):
                    hits += 1
            s2_null1_frac[d].append(hits / N_PERM)

        # S2 null 2: shuffled
        for d in D:
            if d == 0 or s2sz[d] <= s2sz[0]:
                s2_null2_frac[d].append(s2_frac[d][-1])
                continue
            hits = 0
            for _ in range(N_PERM):
                fl2, fl3 = build_shuffled_cascade(
                    sl1, ps_lookup, var_gene_set, None, all_ps_sources, rng)
                pset = (sl1 | fl2) if d == 1 else (sl1 | fl2 | fl3)
                res = run_hypergeom(pset, lib, bg)
                if any(p < SIG_CUTOFF for p in res.values()):
                    hits += 1
            s2_null2_frac[d].append(hits / N_PERM)

        if (idx + 1) % 10 == 0:
            genes = sorted(l1)
            lbl = genes[0] if len(genes) == 1 else "/".join(genes[:2])
            print(f"  [{idx+1}/{n_cres}] {lbl:20s}  "
                  f"tap={real_nsig[0][idx]:2d}/{real_nsig[1][idx]:2d}/"
                  f"{real_nsig[2][idx]:2d}  "
                  f"s2={s2_nsig[0][idx]:2d}/{s2_nsig[1][idx]:2d}/"
                  f"{s2_nsig[2][idx]:2d}")

    # Summarise: % CREs with >= 1 significant annotation
    pct = lambda nsig_list: sum(1 for v in nsig_list if v > 0) / n_cres * 100
    frac_pct = lambda frac_list: np.mean(frac_list) * 100

    # At d=0, nulls didn't run separately — use real terms
    tap_rnd_terms[0]  = tap_ps_terms[0].copy()
    tap_shuf_terms[0] = tap_ps_terms[0].copy()

    return {
        "n_cres": n_cres,
        "pct_tap_ps":   [pct(real_nsig[d]) for d in D],
        "pct_tap_rnd":  [frac_pct(null1_frac[d]) for d in D],
        "pct_tap_shuf": [frac_pct(null2_frac[d]) for d in D],
        "pct_s2_ps":    [frac_pct(s2_frac[d]) for d in D],
        "pct_s2_rnd":   [frac_pct(s2_null1_frac[d]) for d in D],
        "pct_s2_shuf":  [frac_pct(s2_null2_frac[d]) for d in D],
        # Unique term counts
        "uniq_tap_ps":   [len(tap_ps_terms[d]) for d in D],
        "uniq_tap_rnd":  [len(tap_rnd_terms[d]) for d in D],
        "uniq_tap_shuf": [len(tap_shuf_terms[d]) for d in D],
        "uniq_s2_ps":    [len(s2_ps_terms[d]) for d in D],
    }


# ═══════════════════════════════════════════════════════════════════
#  Plot
# ═══════════════════════════════════════════════════════════════════

def plot_panel_A(stats, n_lib_terms):
    print("\n=== Drawing figure ===")

    c_ps, c_rnd, c_shuf = "#E74C3C", "#95A5A6", "#3498DB"
    c_s2 = "#F1948A"

    w = 0.6
    spacing = 1.0
    group_gap = 2.0

    # Group 0: 2 bars (TAP-seq, all predictions)
    g0_xs = [0, spacing]
    # Group 1: 3 bars (real, null1, null2) — TAP-seq only
    g1_xs = [g0_xs[-1] + group_gap + i * spacing for i in range(3)]
    # Group 2: 3 bars (real, null1, null2) — TAP-seq only
    g2_xs = [g1_xs[-1] + group_gap + i * spacing for i in range(3)]

    s = stats

    # ── Values ──
    # Top panel: annotation rate (%)
    pct_g0 = [s["pct_tap_ps"][0], s["pct_s2_ps"][0]]
    pct_g1 = [s["pct_tap_ps"][1], s["pct_tap_rnd"][1], s["pct_tap_shuf"][1]]
    pct_g2 = [s["pct_tap_ps"][2], s["pct_tap_rnd"][2], s["pct_tap_shuf"][2]]
    # Bottom panel: unique signatures
    uniq_g0 = [s["uniq_tap_ps"][0], s["uniq_s2_ps"][0]]
    uniq_g1 = [s["uniq_tap_ps"][1], s["uniq_tap_rnd"][1], s["uniq_tap_shuf"][1]]
    uniq_g2 = [s["uniq_tap_ps"][2], s["uniq_tap_rnd"][2], s["uniq_tap_shuf"][2]]

    # ── Bar styles ──
    g0_cols  = [c_ps, c_s2]
    g0_edges = ["white", "#666"]
    g0_lws   = [0.5, 0.8]
    gh_cols  = [c_ps, c_rnd, c_shuf]
    gh_edges = ["white", "white", "white"]
    gh_lws   = [0.5, 0.5, 0.5]

    g0_lbl_cols = [c_ps, "#B03A2E"]
    gh_lbl_cols = [c_ps, c_rnd, c_shuf]
    gh_bolds    = ["bold", "normal", "normal"]

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, figsize=(10, 8), facecolor="white",
        gridspec_kw={"hspace": 0.35})

    def draw_bars(ax, g0_v, g1_v, g2_v, ylabel, fmt_str, show_xlabels):
        for x, v, c, e, lw in zip(g0_xs, g0_v, g0_cols, g0_edges, g0_lws):
            ax.bar(x, v, w, color=c, edgecolor=e, linewidth=lw, zorder=3)
        for x, v, c, e, lw in zip(g1_xs, g1_v, gh_cols, gh_edges, gh_lws):
            ax.bar(x, v, w, color=c, edgecolor=e, linewidth=lw, zorder=3)
        for x, v, c, e, lw in zip(g2_xs, g2_v, gh_cols, gh_edges, gh_lws):
            ax.bar(x, v, w, color=c, edgecolor=e, linewidth=lw, zorder=3)

        # Value labels
        for x, v, col in zip(g0_xs, g0_v, g0_lbl_cols):
            ax.text(x, v + 0.3, fmt_str.format(v), ha="center", fontsize=10,
                    fontweight="bold", color=col, rotation=90)
        for xs, vals in [(g1_xs, g1_v), (g2_xs, g2_v)]:
            for x, v, col, fw in zip(xs, vals, gh_lbl_cols, gh_bolds):
                ax.text(x, v + 0.3, fmt_str.format(v), ha="center", fontsize=9,
                        fontweight=fw, color=col, rotation=90)

        # X-axis
        if show_xlabels:
            ax.set_xticks([np.mean(g0_xs), np.mean(g1_xs), np.mean(g2_xs)])
            ax.set_xticklabels([
                "CRE \u2192 gene",
                "CRE \u2192 gene \u2192 gene\n(hop 1)",
                "CRE \u2192 gene \u2192 gene \u2192 gene\n(hop 2)",
            ], fontsize=11)
        else:
            ax.set_xticks([np.mean(g0_xs), np.mean(g1_xs), np.mean(g2_xs)])
            ax.set_xticklabels([
                "CRE \u2192 gene",
                "CRE \u2192 gene \u2192 gene\n(hop 1)",
                "CRE \u2192 gene \u2192 gene \u2192 gene\n(hop 2)",
            ], fontsize=11)

        ax.set_ylabel(ylabel, fontsize=13)
        all_v = list(g0_v) + list(g1_v) + list(g2_v)
        ax.set_ylim(0, max(all_v) * 1.35)
        ax.set_xlim(g0_xs[0] - 1, g2_xs[-1] + 1)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", alpha=0.2, zorder=0)

    # ── Top panel: annotation rate ──
    draw_bars(ax_top, pct_g0, pct_g1, pct_g2,
              ylabel="% CREs with \u22651 significant\nCytopus annotation",
              fmt_str="{:.0f}%", show_xlabels=False)

    ax_top.legend(handles=[
        Patch(facecolor=c_ps,   edgecolor="white", label="Real cascade"),
        Patch(facecolor=c_rnd,  edgecolor="white", label="Null 1: random genes"),
        Patch(facecolor=c_shuf, edgecolor="white", label="Null 2: shuffled cascade"),
        Patch(facecolor=c_s2,   edgecolor="#666", linewidth=0.8,
              label="All predictions (S2)"),
    ], fontsize=10.5, loc="upper left", framealpha=0.9)

    # ── Bottom panel: unique signatures ──
    draw_bars(ax_bot, uniq_g0, uniq_g1, uniq_g2,
              ylabel=f"Total unique Cytopus\nsignatures (of {n_lib_terms})",
              fmt_str="{:.0f}", show_xlabels=True)

    ax_bot.axhline(n_lib_terms, color="#999", linestyle="--", linewidth=0.8,
                   zorder=1)
    ax_bot.text(g2_xs[-1] + 0.7, n_lib_terms, f"{n_lib_terms} total",
                fontsize=8, color="#999", va="center")

    for fmt in ["png", "pdf", "svg"]:
        fp = OUTPUT_DIR / f"panel_A_standalone.{fmt}"
        fig.savefig(fp, dpi=300, bbox_inches="tight", facecolor="white")
        print(f"  Saved {fp}")
    plt.close(fig)


# ═══════════════════════════════════════════════════════════════════
#  Entry point
# ═══════════════════════════════════════════════════════════════════

def main():
    enh, ps_lookup, var_gene_set, enh_target_genes, s2_lookup = load_data()
    bg = sorted(var_gene_set)

    lib = load_geneset_library_cached("Cytopus_v1.3")
    n_lib_terms = len(lib)
    print(f"  {n_lib_terms} Cytopus gene sets")

    print("\n=== Enrichment analysis ===")
    stats = run_enrichment(enh, ps_lookup, var_gene_set, enh_target_genes,
                           s2_lookup, lib, bg)
    plot_panel_A(stats, n_lib_terms)
    print("\nDone!")


if __name__ == "__main__":
    main()

# %%
