#%%
"""
Compact TF enrichment figure with motif logos (Figure 4E).

Downloads DoRothEA from OmniPath and JASPAR PFMs at runtime.
Input: input/perturbseq_hits_nolfc.csv (FDR-filtered DESeq2 results, self-hits excluded)

Configuration: DoRothEA A+B+C levels, minDeg >= 5, minDoro >= 50.
Output: output/tf_enrichment.pdf, output/tf_enrichment.png
"""

import csv
import io
import json
import math
import os
import sys
import warnings
from collections import defaultdict
from urllib.request import urlopen, Request
from urllib.error import URLError

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import pandas as pd
import logomaker
from scipy import stats as sp_stats

warnings.filterwarnings('ignore')

# Illustrator-compatible
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['svg.fonttype'] = 'none'
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['font.size'] = 7
plt.rcParams['mathtext.default'] = 'regular'

# === Config ===
CSV_PATH = "/g/stegle/schrod/code/TCell/Figure_5_input/perturbseq_hits_nolfc.csv"
OUTPUT_DIR = "/g/stegle/schrod/code/TCell/Figure_4_output"
OUTPUT = os.path.join(OUTPUT_DIR, 'tf_enrichment')
P_THRESH = 0.1
LFC_THRESH = 0.5
MIN_DEG = 5
MIN_DORO = 50
N_GENES_TESTED = 8001   # genes measured per perturbation (constant across screen)

# Colors
C_CONC = '#18974C'
C_DISC = '#E67E22'
C_UNSN = '#B0B0B0'
C_OBS  = '#2C7FB8'
C_EXP  = '#D0D0CE'

# JASPAR matrix IDs for known TFs
JASPAR_IDS = {
    'STAT2': 'MA1623.2', 'RFX5': 'MA0510.3', 'IRF9': 'MA0653.1',
    'SREBF1': 'MA0829.3', 'SRF': 'MA0083.3', 'MYB': 'MA0100.4',
    'STAT3': 'MA0144.3', 'NFKB2': 'MA0778.2', 'HIF1A': 'MA1106.2',
    'SREBF2': 'MA0828.3', 'STAT5B': 'MA1625.2', 'THAP11': 'MA1573.2',
    'NRF1': 'MA0506.1', 'TFDP1': 'MA1122.2',
}


# ================================================================
# Download DoRothEA from OmniPath
# ================================================================
def download_dorothea():
    """Download DoRothEA interactions from OmniPath API."""
    url = ('https://omnipathdb.org/interactions'
           '?datasets=dorothea&fields=dorothea_level'
           '&dorothea_levels=A,B,C,D&genesymbols=yes')
    print('Downloading DoRothEA from OmniPath ...')
    req = Request(url, headers={'User-Agent': 'Python/figure_script'})
    try:
        with urlopen(req, timeout=120) as resp:
            text = resp.read().decode('utf-8')
    except URLError as e:
        print(f'ERROR: Could not download DoRothEA: {e}', file=sys.stderr)
        sys.exit(1)

    reader = csv.DictReader(io.StringIO(text), delimiter='\t')
    rows = list(reader)
    print(f'  Downloaded {len(rows)} DoRothEA interactions')
    return rows


def parse_dorothea(rows, measured_genes):
    """Parse DoRothEA rows into target sets (A+B+C only),
    restricted to genes in the measured transcriptome."""
    doro_targets = defaultdict(set)
    doro_activated = defaultdict(set)
    doro_repressed = defaultdict(set)

    for row in rows:
        lvl = row.get('dorothea_level', '')
        if 'A' not in lvl and 'B' not in lvl and 'C' not in lvl:
            continue
        tf = row['source_genesymbol']
        tgt = row['target_genesymbol']
        if tgt not in measured_genes:
            continue
        doro_targets[tf].add(tgt)
        cs = row.get('consensus_stimulation', '0')
        ci = row.get('consensus_inhibition', '0')
        cs = cs in ('True', 'true', '1', 1, True)
        ci = ci in ('True', 'true', '1', 1, True)
        if cs and not ci:
            doro_activated[tf].add(tgt)
        elif ci and not cs:
            doro_repressed[tf].add(tgt)

    return doro_targets, doro_activated, doro_repressed


# ================================================================
# Fetch JASPAR PFMs
# ================================================================
def fetch_jaspar_pfm(matrix_id):
    """Fetch a single PFM from JASPAR REST API."""
    url = f'https://jaspar.elixir.no/api/v1/matrix/{matrix_id}/?format=json'
    req = Request(url, headers={'User-Agent': 'Python/figure_script'})
    try:
        with urlopen(req, timeout=30) as resp:
            data = json.loads(resp.read().decode('utf-8'))
        pfm = data.get('pfm', {})
        name = data.get('name', matrix_id)
        return {'jaspar_id': matrix_id, 'name': name, 'pfm': pfm}
    except Exception as e:
        print(f'  Warning: Could not fetch JASPAR {matrix_id}: {e}')
        return None


def fetch_all_jaspar(tf_names):
    """Fetch JASPAR PFMs for all TFs that have known IDs."""
    pfm_data = {}
    needed = [tf for tf in tf_names if tf in JASPAR_IDS]
    if not needed:
        return pfm_data
    print(f'Fetching {len(needed)} JASPAR motifs ...')
    for tf in needed:
        mid = JASPAR_IDS[tf]
        result = fetch_jaspar_pfm(mid)
        if result:
            pfm_data[tf] = result
            print(f'  {tf}: {mid} OK')
    return pfm_data


# ================================================================
# Helpers
# ================================================================
def pfm_to_info_matrix(pfm_dict):
    df = pd.DataFrame(pfm_dict)
    df = df.div(df.sum(axis=1), axis=0)
    info = logomaker.transform_matrix(df, from_type='probability', to_type='information')
    return info


def sig_stars(pval):
    if pval < 0.001:
        return '***'
    elif pval < 0.01:
        return '**'
    elif pval < 0.05:
        return '*'
    return ''


# ================================================================
# Main
# ================================================================
def main():
    # Load Perturb-seq data
    hit_list = pd.read_excel("/g/stegle/schrod/code/TCell/Table S12.xlsx")
    hit_list = hit_list[hit_list["KO"] != hit_list["gene_name"]]  # no self-hits


    all_contrasts = set(hit_list['KO'].unique())
    measured_genes = set(hit_list['gene_name'].unique())
    perturb_data = defaultdict(list)
    for _, row in hit_list.iterrows():
        contrast = row['KO']
        variable = row['gene_name']
        lfc = row['LFC']
        adj_p = row['adj_pval']
        perturb_data[contrast].append((variable, lfc, adj_p))


    N = N_GENES_TESTED
    print(f'Genes tested per perturbation: {N}')
    print(f'Perturbations with hits: {len(all_contrasts)}')

    # Download and parse DoRothEA (targets restricted to measured genes)
    doro_rows = download_dorothea()
    doro_targets, doro_activated, doro_repressed = parse_dorothea(doro_rows, measured_genes)

    # Identify TFs: perturbations where the perturbed gene is a DoRothEA TF
    tf_list = []

    for gene in sorted(all_contrasts):
        targets = doro_targets.get(gene, set())
        if len(targets) < MIN_DORO:
            continue
        if gene not in perturb_data:
            continue
        vdata = perturb_data[gene]

        deg_set = set()
        for var, lfc, adj_p in vdata:
            if abs(lfc) > LFC_THRESH:
                deg_set.add(var)

        if len(deg_set) < MIN_DEG:
            continue

        tf_list.append((gene, vdata, deg_set, targets))

    print(f'TFs passing cutoffs (deg>={MIN_DEG}, doro>={MIN_DORO}, A+B+C): {len(tf_list)}')

    # Fetch JASPAR motifs for passing TFs
    tf_names = [t[0] for t in tf_list]
    pfm_data = fetch_all_jaspar(tf_names)

    # Compute per-TF stats
    all_data = []
    for tf_name, vdata, deg_set, known in tf_list:
        act_targets = doro_activated.get(tf_name, set())
        rep_targets = doro_repressed.get(tf_name, set())

        n_conc = 0
        n_disc = 0
        n_unsn = 0

        for var, lfc, adj_p in vdata:
            if abs(lfc) > LFC_THRESH:
                if var in known:
                    if var in act_targets:
                        if lfc < -LFC_THRESH:
                            n_conc += 1
                        else:
                            n_disc += 1
                    elif var in rep_targets:
                        if lfc > LFC_THRESH:
                            n_conc += 1
                        else:
                            n_disc += 1
                    else:
                        n_unsn += 1

        n_known = len(known)
        overlap = len(deg_set & known)

        # Fisher's exact test:
        # N = total genes tested in this perturbation
        # rows: DEG vs not-DEG; cols: DoRothEA target vs not-target
        a = overlap
        b = len(deg_set) - a
        c = n_known - a
        d = N - a - b - c
        odds, pval = sp_stats.fisher_exact([[a, b], [c, d]], alternative='greater')

        expected_frac = (n_known / N) * 100
        observed_frac = (overlap / len(deg_set)) * 100

        n_dir = n_conc + n_disc
        if n_dir > 0:
            binom_p = sp_stats.binomtest(n_conc, n_dir, 0.5, alternative='greater').pvalue
        else:
            binom_p = 1.0

        all_data.append({
            'name': tf_name, 'n_deg': len(deg_set),
            'expected_frac': expected_frac, 'observed_frac': observed_frac,
            'overlap': overlap, 'n_known': n_known, 'odds': odds, 'pval': pval,
            'n_conc': n_conc, 'n_disc': n_disc, 'n_unsn': n_unsn,
            'n_dir': n_dir, 'binom_p': binom_p,
        })
        print(f'  {tf_name:12s}  DEGs={len(deg_set):4d}  doro={n_known:3d}  '
              f'overlap={overlap:2d}  p={pval:.2e}  '
              f'conc={n_conc} disc={n_disc} unsn={n_unsn}  binom_p={binom_p:.2e}')

    # Sort by Fisher p-value
    all_data.sort(key=lambda d: d['pval'])

    # === Figure ===
    n = len(all_data)
    ROW_H = 0.36
    TOP_PAD = 0.35
    BOT_PAD = 0.45
    RIGHT_PAD = 0.08

    fig_h = ROW_H * n + TOP_PAD + BOT_PAD
    fig_w = 5.5

    fig = plt.figure(figsize=(fig_w, fig_h))

    gs = gridspec.GridSpec(
        1, 3, figure=fig,
        width_ratios=[0.18, 0.42, 0.40],
        left=0.105,
        right=1 - RIGHT_PAD / fig_w,
        top=1 - TOP_PAD / fig_h,
        bottom=BOT_PAD / fig_h,
        wspace=0.06
    )

    y_pos = np.arange(n)
    bar_h = 0.32
    fs_label = 5       # font size for bar labels
    fs_name = 7        # font size for TF names
    all_data_rev = list(reversed(all_data))

    # MOTIF COLUMN
    ax_motif_host = fig.add_subplot(gs[0, 0])
    ax_motif_host.set_xlim(0, 1)
    ax_motif_host.set_ylim(-0.6, n - 0.4)
    ax_motif_host.axis('off')

    for i, d in enumerate(all_data_rev):
        tf_name = d['name']
        y_frac = (y_pos[i] - 0.35 + 0.6) / (n - 0.4 + 0.6)
        h_frac = 0.85 / n
        ax_m = ax_motif_host.inset_axes([0.0, y_frac, 1.0, h_frac])

        if tf_name in pfm_data:
            pfm = pfm_data[tf_name]['pfm']
            info_mat = pfm_to_info_matrix(pfm)
            logomaker.Logo(info_mat, ax=ax_m, color_scheme='classic',
                           show_spines=False, baseline_width=0)
            ax_m.set_ylim(0, 2.3)
        else:
            ax_m.text(0.5, 0.5, '–', transform=ax_m.transAxes,
                      ha='center', va='center', fontsize=5, color='#999')

        ax_m.set_xticks([])
        ax_m.set_yticks([])
        for sp in ax_m.spines.values():
            sp.set_visible(False)

        ax_m.text(-0.06, 0.5, tf_name, transform=ax_m.transAxes,
                  ha='right', va='center', fontsize=fs_name, fontweight='bold',
                  family='Arial')

    # LEFT PANEL: Enrichment
    ax_enr = fig.add_subplot(gs[0, 1])

    max_frac = max(max(d['expected_frac'], d['observed_frac']) for d in all_data) * 1.3
    max_frac = max(math.ceil(max_frac), 3)

    exp_vals = [d['expected_frac'] for d in all_data_rev]
    obs_vals = [d['observed_frac'] for d in all_data_rev]
    obs_colors = [C_OBS if d['pval'] < 0.05 else '#A8C8DC' for d in all_data_rev]

    ax_enr.barh(y_pos + bar_h / 2 + 0.02, exp_vals, height=bar_h,
                color=C_EXP, edgecolor='#999', linewidth=0.4, label='Expected', zorder=2)
    for i in range(n):
        ax_enr.barh(y_pos[i] - bar_h / 2 - 0.02, obs_vals[i], height=bar_h,
                    color=obs_colors[i], edgecolor='#1A5276', linewidth=0.4, zorder=2)
    ax_enr.barh([], [], height=bar_h, color=C_OBS, edgecolor='#1A5276',
                linewidth=0.4, label='Observed')

    for i, d in enumerate(all_data_rev):
        ax_enr.text(d['expected_frac'] + max_frac * 0.012,
                    y_pos[i] + bar_h / 2 + 0.02,
                    '{:.1f}%'.format(d['expected_frac']), va='center', ha='left',
                    fontsize=fs_label, color='#999')
        stars = sig_stars(d['pval'])
        obs_label = '{:.1f}%'.format(d['observed_frac'])
        if stars:
            obs_label += ' ' + stars
        color = '#0A5032' if d['pval'] < 0.05 else '#555'
        ax_enr.text(d['observed_frac'] + max_frac * 0.012,
                    y_pos[i] - bar_h / 2 - 0.02,
                    obs_label, va='center', ha='left',
                    fontsize=fs_label, color=color, fontweight='bold')

    ax_enr.set_xlim(0, max_frac)
    ax_enr.set_ylim(-0.6, n - 0.4)
    ax_enr.set_yticks([])
    ax_enr.set_xlabel('% of DEGs that are\nDoRothEA targets', fontsize=7, labelpad=2)
    ax_enr.set_title('Target enrichment', fontsize=8, fontweight='bold',
                      color='#333', pad=4)
    ax_enr.spines['top'].set_visible(False)
    ax_enr.spines['right'].set_visible(False)
    ax_enr.spines['left'].set_visible(False)
    ax_enr.tick_params(axis='y', length=0)
    ax_enr.tick_params(axis='x', labelsize=6)
    ax_enr.legend(fontsize=6, loc='lower right', framealpha=0.9,
                  handlelength=1, handleheight=0.7)

    for i in range(n - 1):
        ax_enr.axhline(y=i + 0.5, color='#E0E0E0', lw=0.3, zorder=0)

    # RIGHT PANEL: Directional concordance
    ax_dir = fig.add_subplot(gs[0, 2])

    max_hits = max(d['n_conc'] + d['n_disc'] + d['n_unsn'] for d in all_data)
    x_lim = max(max_hits * 1.55, 3)

    conc_vals = [d['n_conc'] for d in all_data_rev]
    disc_vals = [d['n_disc'] for d in all_data_rev]
    unsn_vals = [d['n_unsn'] for d in all_data_rev]

    ax_dir.barh(y_pos, conc_vals, height=bar_h * 1.4,
                color=C_CONC, edgecolor='#333', linewidth=0.3,
                label='Concordant', zorder=2)
    ax_dir.barh(y_pos, disc_vals, left=conc_vals, height=bar_h * 1.4,
                color=C_DISC, edgecolor='#333', linewidth=0.3,
                label='Discordant', zorder=2)
    lefts_unsn = [c + d for c, d in zip(conc_vals, disc_vals)]
    ax_dir.barh(y_pos, unsn_vals, left=lefts_unsn, height=bar_h * 1.4,
                color=C_UNSN, edgecolor='#333', linewidth=0.3,
                label='Unsigned', zorder=2)

    for i, d in enumerate(all_data_rev):
        total = d['n_conc'] + d['n_disc'] + d['n_unsn']
        x_text = total + x_lim * 0.02
        stars = sig_stars(d['binom_p'])
        if d['n_dir'] > 0:
            label = '{}/{}'.format(d['n_conc'], d['n_dir'])
            if stars:
                label += ' ' + stars
            color = '#0A5032' if d['binom_p'] < 0.05 else '#555'
        elif total > 0:
            label = '(no dir.)'
            color = '#999'
        else:
            label = '—'
            color = '#999'
        ax_dir.text(x_text, y_pos[i], label,
                    va='center', ha='left', fontsize=fs_label,
                    color=color, fontweight='bold')

    ax_dir.set_xlim(0, x_lim)
    ax_dir.set_ylim(-0.6, n - 0.4)
    ax_dir.set_yticks([])
    ax_dir.set_xlabel('DoRothEA targets\namong DEGs (count)', fontsize=7, labelpad=2)
    ax_dir.set_title('Directional concordance', fontsize=8, fontweight='bold',
                      color='#333', pad=4)
    ax_dir.spines['top'].set_visible(False)
    ax_dir.spines['right'].set_visible(False)
    ax_dir.spines['left'].set_visible(False)
    ax_dir.tick_params(axis='y', length=0)
    ax_dir.tick_params(axis='x', labelsize=6)
    ax_dir.legend(fontsize=6, loc='lower right', framealpha=0.9,
                  handlelength=1, handleheight=0.7)

    for i in range(n - 1):
        ax_dir.axhline(y=i + 0.5, color='#E0E0E0', lw=0.3, zorder=0)

    # Save
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    fig.savefig(OUTPUT + '.pdf', format='pdf', dpi=300, bbox_inches='tight')
    fig.savefig(OUTPUT + '.png', format='png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'\nSaved: {OUTPUT}.pdf')
    print(f'Saved: {OUTPUT}.png')


if __name__ == '__main__':
    main()

# %%
