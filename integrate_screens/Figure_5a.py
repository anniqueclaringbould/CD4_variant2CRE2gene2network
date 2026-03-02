#%%
#!/usr/bin/env python3
"""
Figure 5A — Genome-wide regulatory network sparsity visualization.

Four-layer network ordered by chromosomal position (chr1-22, chrX):
  L0  Enhancer-like CREs (colored diamonds, one per disease)
  L1  TAP-seq tested genes (gray = untargeted, blue = has enhancer hit)
  L2  Variable genes — Perturb-seq hop 1 (gray = inactive, green = active)
  L3  Variable genes — Perturb-seq hop 2 (gray = inactive, gold = active)

Edges drawn as straight lines via LineCollection:
  Red   CRE -> gene         (TAP-seq validated)
  Green gene -> hop 1 gene  (Perturb-seq)
  Gold  hop 1 -> hop 2 gene (Perturb-seq)

Input files (in input/):
  gene_set.csv                                — variably expressed genes
  cre_perturbation_screen_results_traits.csv  — TAP-seq CRE perturbation screen
  perturbseq_hits_nolfc.csv                   — Perturb-seq DESeq2 results
  ncbiRefSeq_cache.txt.gz                     — UCSC RefSeq hg38 coordinates

Output (in output/):
  figure_5A_network.pdf / .svg / .png
"""

import re
import gzip
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
from matplotlib.gridspec import GridSpec

# ── Paths ────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).resolve().parent
INPUT_DIR  = Path("/g/stegle/schrod/code/TCell/Figure_5_input")
OUTPUT_DIR = Path("/g/stegle/schrod/code/TCell/Figure_5_output")

# ── Filtering thresholds ─────────────────────────────────────────────────────
TOP_N       = 100000           # max Perturb-seq targets per source gene
LFC_CUTOFF  = 0.2          # min |log2 fold-change| for Perturb-seq edges
DIST_CUTOFF = 1_000_000    # max |dist_to_tss| in bp for TAP-seq edges

# ── Figure parameters ────────────────────────────────────────────────────────
FIG_W, FIG_H = 7.0, 5.0
Y_L0, Y_L1, Y_L2, Y_L3 = 0.90, 0.66, 0.38, 0.10   # layer y-positions
BAND_H  = 0.055    # half-height of each tick band
CHR_GAP = 0.003    # normalized gap between chromosome columns

COLOR = {0: "#C0392B", 1: "#2471A3", 2: "#27AE60", 3: "#D4AC0D"}
MULTI_DISEASE = "#2C3E50"

CHROM_ORDER = [f"chr{i}" for i in range(1, 23)] + ["chrX"]
CHROM_SET   = set(CHROM_ORDER)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})


# ═════════════════════════════════════════════════════════════════════════════
# 1. LOAD & FILTER DATA
# ═════════════════════════════════════════════════════════════════════════════

print("Loading variable gene list...")
var_genes = sorted(
    pd.read_csv(INPUT_DIR / "gene_set.csv")["bpcells_name"].unique()
)
var_gene_set = set(var_genes)
n_var = len(var_genes)
print(f"  {n_var} variable genes")

# ── TAP-seq ──────────────────────────────────────────────────────────────────
print("Loading TAP-seq data...")
tapseq = pd.read_csv(
    INPUT_DIR / "cre_perturbation_screen_results_traits.csv", low_memory=False
)

enh = tapseq[
    (tapseq["enhancer_like_interaction"] == True)
    & (tapseq["significant"] == True)
].copy()
enh = enh[enh["dist_to_tss"].abs() <= DIST_CUTOFF]
print(f"  Enhancer-like edges after dist filter: {len(enh)}")

# All genes tested across any panel (for L1 background ticks)
gene_panel = tapseq.groupby("response_id")["panel"].min().to_dict()
all_tested_genes = sorted(gene_panel, key=lambda g: (gene_panel[g], g))
enh_target_genes = set(enh["response_id"].unique())
print(f"  Tested genes: {len(all_tested_genes)}, with enhancer hits: {len(enh_target_genes)}")

ps = pd.read_excel("/g/stegle/schrod/code/TCell/Table S12.xlsx")
ps.rename(columns={"KO": "contrast", "gene_name": "variable", "adj_pval": "adj_p_value", "LFC": "log_fc"}, inplace=True)
ps = ps[ps["contrast"] != ps["variable"]]                          # no self
ps = ps[ps["log_fc"].abs() >= LFC_CUTOFF]                          # LFC gate
ps = ps[ps["contrast"].isin(var_gene_set) & ps["variable"].isin(var_gene_set)]  # both var
ps = ps.sort_values("adj_p_value")
ps["rank"] = ps.groupby("contrast").cumcount() + 1
ps = ps[ps["rank"] <= TOP_N]                                       # top-N cap
print(f"  Filtered Perturb-seq edges: {len(ps):,}")

# Build lookup: source gene -> list of targets
ps_lookup = {}
for _, r in ps.iterrows():
    ps_lookup.setdefault(r["contrast"], []).append(r["variable"])


# ═════════════════════════════════════════════════════════════════════════════
# 2. BUILD NETWORK EDGES (3 layers)
# ═════════════════════════════════════════════════════════════════════════════

print("Building network edges...")

# L0 -> L1  (TAP-seq: CRE -> gene)
tapseq_edges = list(zip(enh["grna_target"], enh["response_id"]))

# L1 -> L2  (Perturb-seq hop 1)
gene_active_l2 = set()
l1_to_l2_edges = []
for gene in enh_target_genes:
    for tgt in ps_lookup.get(gene, []):
        l1_to_l2_edges.append((gene, tgt))
        gene_active_l2.add(tgt)

# L2 -> L3  (Perturb-seq hop 2)
gene_active_l3 = set()
l2_to_l3_edges = []
for gene in gene_active_l2:
    for tgt in ps_lookup.get(gene, []):
        l2_to_l3_edges.append((gene, tgt))
        gene_active_l3.add(tgt)

print(f"  L0->L1: {len(tapseq_edges)}  L1->L2: {len(l1_to_l2_edges)}  L2->L3: {len(l2_to_l3_edges)}")


# ═════════════════════════════════════════════════════════════════════════════
# 3. DISEASE COLOR MAPPING FOR CREs
# ═════════════════════════════════════════════════════════════════════════════

def _clean(d):
    d = d.strip().replace("_", " ")
    return "Primary biliary cholangitis" if d == "Primary biliary cirrhosis" else d

cre_diseases = {}
for _, r in enh.iterrows():
    cre = r["grna_target"]
    cre_diseases.setdefault(cre, set()).update(
        _clean(d) for d in str(r["trait"]).split(",")
    )

all_diseases = sorted({d for ds in cre_diseases.values() for d in ds})
cmap = plt.cm.tab20
disease_colors = {d: cmap(i / len(all_diseases)) for i, d in enumerate(all_diseases)}


# ═════════════════════════════════════════════════════════════════════════════
# 4. RESOLVE GENOMIC COORDINATES
# ═════════════════════════════════════════════════════════════════════════════

print("Resolving genomic coordinates...")

# 4a. Gene coords from TAP-seq columns
gene_chrom = {}
for _, r in tapseq.drop_duplicates("response_id").iterrows():
    g, c = r["response_id"], str(r.get("gene_chr", ""))
    if c in CHROM_SET:
        gene_chrom[g] = (c, int(r["gene_start"]))

# 4b. RefSeq fallback
print("  Loading RefSeq cache...")
with gzip.open(INPUT_DIR / "ncbiRefSeq_cache.txt.gz", "rt") as f:
    refseq_data = f.read()

refseq = {}
for line in refseq_data.strip().split("\n"):
    cols = line.split("\t")
    if len(cols) < 13:
        continue
    chrom, start, end, sym = cols[2], int(cols[4]), int(cols[5]), cols[12]
    if chrom not in CHROM_SET:
        continue
    if sym not in refseq or (end - start) > (refseq[sym][1] - refseq[sym][0]):
        refseq[sym] = (start, end, chrom)

for g in set(var_genes) | set(all_tested_genes):
    if g not in gene_chrom and g in refseq:
        s, e, c = refseq[g]
        gene_chrom[g] = (c, s)

mapped = sum(1 for g in var_genes if g in gene_chrom)
print(f"  Variable genes mapped: {mapped}/{n_var}")

# 4c. CRE coords (parse name, fallback to TAP-seq columns)
all_cres = sorted(enh["grna_target"].unique())
cre_chrom = {}
for cre in all_cres:
    m = re.match(r"(chr\w+):(\d+)-(\d+)", cre)
    if m and m.group(1) in CHROM_SET:
        cre_chrom[cre] = (m.group(1), (int(m.group(2)) + int(m.group(3))) // 2)

for _, r in enh.drop_duplicates("grna_target").iterrows():
    cre = r["grna_target"]
    if cre not in cre_chrom:
        c = str(r.get("pert_chr", ""))
        if c in CHROM_SET:
            cre_chrom[cre] = (c, int(r["pert_start"]))

print(f"  CREs mapped: {len(cre_chrom)}/{len(all_cres)}")


# ═════════════════════════════════════════════════════════════════════════════
# 5. CHROMOSOME-PROPORTIONAL LAYOUT
# ═════════════════════════════════════════════════════════════════════════════

print("Computing chromosome layout...")

def _assign_x(items, chr_starts, chr_widths):
    """Assign x-positions within chromosome columns for (name, chrom, pos) items."""
    by_chr = {c: [] for c in CHROM_ORDER}
    for name, chrom, pos in items:
        if chrom in by_chr:
            by_chr[chrom].append((pos, name))
    for c in CHROM_ORDER:
        by_chr[c].sort()

    xmap = {}
    for c in CHROM_ORDER:
        genes = by_chr[c]
        n = len(genes)
        if n == 0:
            continue
        xs, w = chr_starts[c], chr_widths[c]
        for i, (_, name) in enumerate(genes):
            xmap[name] = xs + w / 2 if n == 1 else xs + w * i / (n - 1)
    return xmap

# Column widths from variable genes (L2/L3 define the layout)
var_mapped = [(g, *gene_chrom[g]) for g in var_genes if g in gene_chrom]
chr_counts = {c: 0 for c in CHROM_ORDER}
for _, chrom, _ in var_mapped:
    chr_counts[chrom] += 1

total = sum(chr_counts.values())
usable = 1.0 - CHR_GAP * (len(CHROM_ORDER) - 1)

chr_widths, chr_starts = {}, {}
cursor = 0.0
for c in CHROM_ORDER:
    w = usable * chr_counts[c] / total if total > 0 else 0.0
    chr_widths[c] = w
    chr_starts[c] = cursor
    cursor += w + CHR_GAP

# Assign x-positions per layer
var_gene_x = _assign_x(var_mapped, chr_starts, chr_widths)

tested_mapped = [(g, *gene_chrom[g]) for g in all_tested_genes if g in gene_chrom]
tested_gene_x = _assign_x(tested_mapped, chr_starts, chr_widths)

cre_mapped = [(cre, *cre_chrom[cre]) for cre in all_cres if cre in cre_chrom]
cre_x = _assign_x(cre_mapped, chr_starts, chr_widths)

n_var_mapped = len(var_gene_x)
print(f"  CREs: {len(cre_x)}, tested genes: {len(tested_gene_x)}, var genes: {n_var_mapped}")


# ═════════════════════════════════════════════════════════════════════════════
# 6. DRAW FIGURE
# ═════════════════════════════════════════════════════════════════════════════

print("Drawing figure...")

fig = plt.figure(figsize=(FIG_W, FIG_H), dpi=300, facecolor="white")
gs = GridSpec(2, 1, figure=fig, height_ratios=[0.82, 0.18], hspace=0.08)
ax     = fig.add_subplot(gs[0])
ax_leg = fig.add_subplot(gs[1])
ax.set_facecolor("white")
ax_leg.set_facecolor("white")
ax_leg.axis("off")

layer_ys = [Y_L0, Y_L1, Y_L2, Y_L3]

# ── Chromosome background shading & separators ──────────────────────────────
for i, c in enumerate(CHROM_ORDER):
    if chr_counts[c] == 0:
        continue
    xl = chr_starts[c] - CHR_GAP * 0.3
    xr = chr_starts[c] + chr_widths[c] + CHR_GAP * 0.3
    if i % 2 == 0:
        for y in layer_ys:
            ax.fill_between([xl, xr], y - BAND_H, y + BAND_H,
                            color="#F0F1F2", zorder=0, alpha=0.7)
    if i > 0:
        xsep = chr_starts[c] - CHR_GAP / 2
        for y in layer_ys:
            ax.plot([xsep, xsep], [y - BAND_H, y + BAND_H],
                    color="#D5D8DC", lw=0.3, zorder=1, alpha=0.6)

# ── Chromosome labels ────────────────────────────────────────────────────────
for c in CHROM_ORDER:
    if chr_counts[c] == 0:
        continue
    cx = chr_starts[c] + chr_widths[c] / 2
    label = c.replace("chr", "")
    w = chr_widths[c]
    fs, rot = (4.5, 0) if w > 0.025 else (3.5, 0) if w > 0.012 else (3.5, 90)
    ax.text(cx, 1.02, label, fontsize=fs, ha="center", va="bottom",
            color="#2C3E50", fontweight="bold", rotation=rot,
            transform=ax.get_xaxis_transform(), clip_on=False)

ax.axhline(y=0.985, color="#2C3E50", lw=0.6, zorder=2)
ax.text(-0.008, 0.985, "Chr", fontsize=4.5, ha="right", va="center",
        color="#2C3E50", fontweight="bold")

# ── Layer labels ─────────────────────────────────────────────────────────────
lx = -0.01
n_l1_active = len(enh_target_genes)
n_l2_active = len(gene_active_l2)
n_l3_active = len(gene_active_l3)

ax.text(lx, Y_L0, f"CREs\n(n={len(cre_x)})", fontsize=5, fontweight="bold",
        color=COLOR[0], ha="right", va="center", linespacing=1.1)
ax.text(lx, Y_L1, f"Tested genes\n({n_l1_active}/{len(tested_gene_x)})",
        fontsize=5, fontweight="bold", color=COLOR[1], ha="right", va="center",
        linespacing=1.1)
ax.text(lx, Y_L2, f"Hop 1\n({n_l2_active}/{n_var_mapped})", fontsize=5,
        fontweight="bold", color=COLOR[2], ha="right", va="center", linespacing=1.1)
ax.text(lx, Y_L3, f"Hop 2\n({n_l3_active}/{n_var_mapped})", fontsize=5,
        fontweight="bold", color=COLOR[3], ha="right", va="center", linespacing=1.1)

# ── Edges (bottom-up: gold -> green -> red) ──────────────────────────────────
def _edge_segs(edges, src_x, tgt_x, y_src, y_tgt):
    return [(( src_x[s], y_src), (tgt_x[t], y_tgt))
            for s, t in edges if s in src_x and t in tgt_x]

def _edge_segs_with_limit(edges, src_x, tgt_x, y_src, y_tgt, max_visible=20):
    """Return segments split into visible (≤20 per source) and invisible (>20)."""
    visible, invisible = [], []
    src_current = {}
    sources_with_connections = set()
    
    for s, t in edges:
        if s not in src_x or t not in tgt_x:
            continue
        
        sources_with_connections.add(s)
        src_current[s] = src_current.get(s, 0) + 1
        seg = ((src_x[s], y_src), (tgt_x[t], y_tgt))
        
        if src_current[s] <= max_visible:
            visible.append(seg)
        else:
            invisible.append(seg)
    
    return visible, invisible, len(sources_with_connections)

segs_vis_l2l3, segs_invis_l2l3, n_src_l2l3 = _edge_segs_with_limit(l2_to_l3_edges, var_gene_x, var_gene_x, Y_L2, Y_L3)
if segs_vis_l2l3:
    ax.add_collection(LineCollection(segs_vis_l2l3, colors=COLOR[3], linewidths=0.2, alpha=0.08, zorder=2))
if segs_invis_l2l3:
    ax.add_collection(LineCollection(segs_invis_l2l3, colors=COLOR[3], linewidths=0.0, alpha=0.08, zorder=2))

segs_vis_l1l2, segs_invis_l1l2, n_src_l1l2 = _edge_segs_with_limit(l1_to_l2_edges, tested_gene_x, var_gene_x, Y_L1, Y_L2)
if segs_vis_l1l2:
    ax.add_collection(LineCollection(segs_vis_l1l2, colors=COLOR[2], linewidths=0.3, alpha=0.12, zorder=3))
if segs_invis_l1l2:
    ax.add_collection(LineCollection(segs_invis_l1l2, colors=COLOR[2], linewidths=0.0, alpha=0.12, zorder=3))

segs = _edge_segs(tapseq_edges, cre_x, tested_gene_x, Y_L0, Y_L1)
n_src_l0l1 = len(set(s for s, t in tapseq_edges if s in cre_x and t in tested_gene_x))
if segs:
    ax.add_collection(LineCollection(segs, colors=COLOR[0], linewidths=0.5, alpha=0.30, zorder=4))

# Print edge statistics
total_l1l2 = len(segs_vis_l1l2) + len(segs_invis_l1l2)
total_l2l3 = len(segs_vis_l2l3) + len(segs_invis_l2l3)
total_l0l1 = len(segs)

# Calculate how many CRE-targeted genes have Perturb-seq edges
cre_targeted_with_perturbseq = len(set(s for s, t in l1_to_l2_edges))

print(f"\nEdge statistics:")
print(f"  L0→L1 (TAP-seq):    {total_l0l1:,} edges shown (no limit), {n_src_l0l1} source CREs")
print(f"  L1→L2 (hop 1):      {len(segs_vis_l1l2):,} / {total_l1l2:,} edges shown ({len(segs_invis_l1l2):,} hidden), {n_src_l1l2} source genes")
print(f"                      ({cre_targeted_with_perturbseq}/{len(enh_target_genes)} CRE-targeted genes have Perturb-seq edges)")
print(f"  L2→L3 (hop 2):      {len(segs_vis_l2l3):,} / {total_l2l3:,} edges shown ({len(segs_invis_l2l3):,} hidden), {n_src_l2l3} source genes")
print(f"  Total visible:      {total_l0l1 + len(segs_vis_l1l2) + len(segs_vis_l2l3):,} / {total_l0l1 + total_l1l2 + total_l2l3:,}")

# ── Ticks (inactive gray + active colored) ──────────────────────────────────
def _draw_ticks(ax, xmap, y, active, color, z):
    inactive, active_segs = [], []
    for g, x in xmap.items():
        seg = [(x, y - BAND_H), (x, y + BAND_H)]
        (active_segs if g in active else inactive).append(seg)
    if inactive:
        ax.add_collection(LineCollection(inactive, colors="#D5D8DC", linewidths=0.2,
                                         alpha=0.5, zorder=z))
    if active_segs:
        ax.add_collection(LineCollection(active_segs, colors=color, linewidths=0.6,
                                         alpha=0.90, zorder=z + 1))

_draw_ticks(ax, var_gene_x,    Y_L3, gene_active_l3,  COLOR[3], z=5)
_draw_ticks(ax, var_gene_x,    Y_L2, gene_active_l2,  COLOR[2], z=8)
_draw_ticks(ax, tested_gene_x, Y_L1, enh_target_genes, COLOR[1], z=11)

# ── CRE diamonds ─────────────────────────────────────────────────────────────
cre_xs, cre_cols = [], []
for cre in [c for c, _, _ in cre_mapped]:
    if cre in cre_x:
        cre_xs.append(cre_x[cre])
        ds = cre_diseases.get(cre, set())
        cre_cols.append(disease_colors[list(ds)[0]] if len(ds) == 1 else MULTI_DISEASE)

ax.scatter(cre_xs, [Y_L0] * len(cre_xs), s=12, c=cre_cols, marker="D",
           alpha=0.95, linewidths=0.3, edgecolors="white", zorder=14)

# ── Layer separators ─────────────────────────────────────────────────────────
for ys in [(Y_L0 + Y_L1) / 2, (Y_L1 + Y_L2) / 2, (Y_L2 + Y_L3) / 2]:
    ax.axhline(y=ys, color="#ECF0F1", lw=0.4, ls="--", zorder=0)

# ── Axes ─────────────────────────────────────────────────────────────────────
ax.set_xlim(-0.008, 1.008)
ax.set_ylim(-0.03, 1.02)
ax.set_xticks([])
ax.set_yticks([])
for sp in ax.spines.values():
    sp.set_visible(False)


# ═════════════════════════════════════════════════════════════════════════════
# 7. LEGEND
# ═════════════════════════════════════════════════════════════════════════════

handles = []

# Disease colors
for d in all_diseases:
    handles.append(Line2D([0], [0], marker="D", color="white",
                          markerfacecolor=disease_colors[d], markersize=3.5,
                          linewidth=0, label=d))
handles.append(Line2D([0], [0], marker="D", color="white",
                       markerfacecolor=MULTI_DISEASE, markersize=3.5,
                       linewidth=0, label="Multiple diseases"))

handles.append(Line2D([0], [0], color="white", linewidth=0, label=" "))

# Edge types
handles.append(Line2D([0], [0], color=COLOR[0], lw=1.5, alpha=0.5,
                       label=f"CRE \u2192 gene (TAP-seq, n={len(tapseq_edges)})"))
handles.append(Line2D([0], [0], color=COLOR[2], lw=1.2, alpha=0.4,
                       label=f"Gene \u2192 hop 1 (PS, n={len(l1_to_l2_edges)})"))
handles.append(Line2D([0], [0], color=COLOR[3], lw=1.2, alpha=0.4,
                       label=f"Hop 1 \u2192 hop 2 (PS, n={len(l2_to_l3_edges)})"))

handles.append(Line2D([0], [0], color="white", linewidth=0, label=" "))

# Node states
handles.append(Line2D([0], [0], marker="|", color=COLOR[2], markersize=5,
                       markeredgewidth=1.2, linewidth=0, label="Active (has connections)"))
handles.append(Line2D([0], [0], marker="|", color="#D5D8DC", markersize=4,
                       markeredgewidth=0.8, linewidth=0, label="Inactive"))

leg = ax_leg.legend(
    handles=handles, loc="center", fontsize=4.5, ncol=4,
    frameon=True, fancybox=True, framealpha=0.95, edgecolor="#BDC3C7",
    handletextpad=0.4, columnspacing=1.0, borderpad=0.8, labelspacing=0.35,
    title=(f"CRE color = disease association | Thresholds: "
           f"adj_p < 0.1, |LFC| \u2265 {LFC_CUTOFF}, top {TOP_N}/gene"),
    title_fontsize=5,
)
leg._legend_box.align = "left"


# ═════════════════════════════════════════════════════════════════════════════
# 8. TITLE & SAVE
# ═════════════════════════════════════════════════════════════════════════════

ax.set_title(
    "Regulatory network sparsity\n"
    f"CRE \u2192 gene \u2192 variable gene cascades | "
    f"{n_var_mapped:,} genes ordered by genomic position",
    fontsize=7, fontweight="bold", pad=12, color="#2C3E50", linespacing=1.3,
)

plt.subplots_adjust(left=0.09, right=0.995, top=0.89, bottom=0.02)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
for ext, dpi in [("pdf", 300), ("svg", None), ("png", 600)]:
    p = OUTPUT_DIR / f"figure_5A_network.{ext}"
    kw = {"bbox_inches": "tight", "facecolor": "white"}
    if dpi:
        kw["dpi"] = dpi
    plt.savefig(p, **kw)
    print(f"  Saved: {p}")

plt.close()
print("Done.")


# %%
