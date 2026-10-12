"""Publication figures (PNG + PDF + SVG, English) for the reptile and amphibian study."""
import json

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import RESULTS, ROOT

OUT = ROOT / 'paper' / 'figures'
OUT.mkdir(parents=True, exist_ok=True)
INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
GROUPS = {'Testudines': ('#2a78d6', 'Turtles'), 'Crocodylia': ('#eb6834', 'Crocodilians'),
          'Lepidosauria': ('#1baf7a', 'Lizards and snakes'), 'Amphibia': ('#4a3aa7', 'Amphibians')}
BLUE, ORANGE = '#2a78d6', '#eb6834'
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Segoe UI', 'DejaVu Sans'], 'font.size': 9,
                     'axes.edgecolor': AXIS, 'axes.labelcolor': INK2, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.titlesize': 10, 'axes.titleweight': 'bold', 'axes.titlecolor': INK,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': 0.6, 'legend.frameon': False, 'figure.dpi': 150, 'savefig.dpi': 300,
                     'savefig.bbox': 'tight', 'pdf.fonttype': 42, 'svg.fonttype': 'none'})


def save(fig, name):
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(OUT / f'{name}.{ext}')
    plt.close(fig)


def by_group(ax, x, y, st, s=30):
    for g, (col, lab) in GROUPS.items():
        m = (st.group == g).to_numpy()
        if m.any():
            ax.scatter(x[m], y[m], s=s, color=col, edgecolor='white', linewidth=1.1, zorder=3, label=f'{lab} (n={m.sum()})')


def label_species(ax, x, y, st, names):
    for k, off in names.items():
        i = np.where(st.binomial.to_numpy() == k)[0]
        if len(i):
            ax.annotate(k.replace('_', ' ').capitalize(), (x[i[0]], y[i[0]]), xytext=off, textcoords='offset points',
                        fontsize=7, color=INK2, style='italic')


def residual_longevity(st):
    x, y = np.log10(st.mass_g.to_numpy()), np.log10(st.max_lifespan_yr.to_numpy())
    b = np.polyfit(x, y, 1)
    return y - np.polyval(b, x)


def fig_species(st):
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    x, y = st.mass_g.to_numpy(), st.max_lifespan_yr.to_numpy()
    by_group(ax, x, y, st, 34)
    label_species(ax, x, y, st, {'emys_orbicularis': (5, 3), 'anolis_carolinensis': (5, -9), 'xenopus_laevis': (5, 3),
                                 'alligator_mississippiensis': (5, 3), 'chelonia_mydas': (5, -9)})
    ax.set(xscale='log', yscale='log', xlabel='Adult body mass (g)', ylabel='Maximum recorded lifespan (years)',
           title=f'{len(st)} reptiles and amphibians analysed')
    ax.legend(loc='upper left', fontsize=7.5)
    save(fig, 'fig1_species')


def fig_controls(c):
    fig, axes = plt.subplots(1, 3, figsize=(8.4, 2.9))
    axes[0].barh([1, 0], [c['C1_mean_observed_ddG'], c['C1_mean_alternative_ddG']], color=[BLUE, AXIS], height=0.5)
    axes[0].set_yticks([1, 0], ['Observed\ndifferences', 'Possible\nalternatives'])
    axes[0].set(xlabel='Mean predicted ΔΔG (kcal/mol)', title='C1 · purifying selection')
    axes[1].barh([1, 0], [c['C2_mean_abs_alt_buried'], c['C2_mean_abs_alt_exposed']], color=[BLUE, AXIS], height=0.5)
    axes[1].set_yticks([1, 0], ['Buried\nsites', 'Exposed\nsites'])
    axes[1].set(xlabel='Mean |ΔΔG| of alternatives', title='C2 · burial')
    vals = [c['C3_helix_region_to_P'], c['C3_helix_region_to_G'], c['C3_helix_region_mean_alt'],
            c['C3_strand_region_to_P'], c['C3_strand_region_to_G'], c['C3_strand_region_mean_alt']]
    axes[2].barh(range(6)[::-1], vals, height=0.6, color=[BLUE, BLUE, AXIS, BLUE, BLUE, AXIS])
    axes[2].set_yticks(range(6)[::-1], ['Helix → Pro', 'Helix → Gly', 'Helix (mean)', 'Strand → Pro', 'Strand → Gly', 'Strand (mean)'])
    axes[2].set(xlabel='Mean predicted ΔΔG', title='C3 · Pro/Gly in structure')
    for a in axes:
        a.grid(axis='y', visible=False)
    fig.tight_layout()
    save(fig, 'fig2_controls')


def fig_primary(st, s):
    p = s['P1']
    x = ((st.M_all - st.M_all.mean()) / st.M_all.std(ddof=1)).to_numpy()
    y = residual_longevity(st)
    fig, ax = plt.subplots(figsize=(5.6, 4.1))
    by_group(ax, x, y, st)
    label_species(ax, x, y, st, {'emys_orbicularis': (5, 3), 'anolis_carolinensis': (5, -9), 'pyxicephalus_adspersus': (5, 3)})
    ax.set(xlabel='Mean ΔΔG per difference (SD; all genes; zebrafish template)',
           ylabel='Lifespan relative to mass (log10 residual)', title='P1 · predicted stability and lifespan')
    ax.text(0.02, 0.03, f"PGLS: β = {p['beta']:.3f} log10 years/SD\n95% CI {p['ci95_low']:.3f} to {p['ci95_high']:.3f}\n"
            f"bootstrap P = {p['p_boot']:.3f} · λ = {p['lambda']:.2f} · n = {p['n']}", transform=ax.transAxes,
            ha='left', va='bottom', fontsize=7.5, color=INK2)
    ax.legend(loc='upper right', fontsize=7.5)
    save(fig, 'fig3_primary')


def fig_competitive(s):
    c = s['P2_competitive']
    null = np.load(RESULTS / 'P2_null_betas.npy')
    fig, ax = plt.subplots(figsize=(5.4, 3.2))
    ax.hist(null, bins=40, color=AXIS, edgecolor='white', linewidth=0.8)
    ax.axvline(c['beta_genage'], color=ORANGE, lw=2)
    ax.annotate(f"GenAge (n={c['set_size']})\nempirical P = {c['p_empirical_two_sided']:.3f}",
                (c['beta_genage'], ax.get_ylim()[1] * 0.9), xytext=(6, 0), textcoords='offset points', fontsize=7.5, color=INK2)
    ax.set(xlabel='β of ΔΔG per difference (log10 years/SD)', ylabel='Random background sets',
           title=f"P2 · GenAge vs {c['n_null_sets']} background sets")
    ax.grid(axis='x', visible=False)
    save(fig, 'fig4_competitive')


def fig_sensitivity(s):
    sens = pd.read_csv(RESULTS / 'sensitivity_P1.csv')
    sens = sens[sens.beta.notna()]
    names = {'lambda_fixed_1': 'λ fixed at 1', 'variant_no_mask': 'No window mask', 'variant_plddt90': 'pLDDT ≥ 90 only',
             'variant_buried': 'Buried sites only', 'variant_exposed': 'Exposed sites only', 'reptiles_only': 'Reptiles only',
             'amphibians_only': 'Amphibians only', 'without_turtles': 'Without turtles', 'squamates_only': 'Squamates only',
             'anage_mass_only': 'AnAge masses only', 'no_mass_covariate': 'No mass covariate',
             'adjust_log10_R': 'Adjusted for divergence from ancestor', 'background_genes_only': 'Background genes only',
             'reptiles_adjust_body_temp_c': 'Reptiles, + body temperature',
             'reptiles_adjust_range_mean_temp_c': 'Reptiles, + range temperature', 'OLS_HC3': 'OLS, HC3 errors'}
    rows = [('Primary (P1)', s['P1']['beta'], s['P1']['ci95_low'], s['P1']['ci95_high'])]
    rows += [(names.get(r.analysis, r.analysis), r.beta, r.ci95_low, r.ci95_high) for r in sens.itertuples()]
    fig, ax = plt.subplots(figsize=(5.8, 0.28 * len(rows) + 0.8))
    ys = np.arange(len(rows))[::-1]
    for k, (yy, (n, b, lo, hi)) in enumerate(zip(ys, rows)):
        col = ORANGE if k == 0 else BLUE
        ax.plot([lo, hi], [yy, yy], color=col, lw=2, solid_capstyle='round')
        ax.scatter([b], [yy], color=col, s=26, zorder=3, edgecolor='white', linewidth=1)
    ax.axvline(0, color=MUTED, lw=1, ls='--')
    ax.set_yticks(ys, [r[0] for r in rows])
    ax.set(xlabel='β (log10 years per SD) and approximate 95% CI', title='Sensitivity of P1')
    ax.grid(axis='y', visible=False)
    save(fig, 'fig5_sensitivity')


def fig_genes():
    g = pd.read_csv(RESULTS / 'per_gene_exploratory.csv')
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    bg, ga = g[g.gene_set == 'background'], g[g.gene_set == 'genage']
    ax.scatter(bg.beta, -np.log10(bg.p_wald), s=9, color=AXIS, label=f'Background (n={len(bg)})', zorder=2)
    ax.scatter(ga.beta, -np.log10(ga.p_wald), s=14, color=ORANGE, label=f'GenAge (n={len(ga)})', zorder=3)
    if (g.q_bh <= 0.05).any():
        thr = g.loc[g.q_bh <= 0.05, 'p_wald'].max()
        ax.axhline(-np.log10(thr), color=MUTED, lw=1, ls='--')
        ax.text(ax.get_xlim()[0], -np.log10(thr), ' FDR 5%', va='bottom', fontsize=7, color=MUTED)
    ax.set(xlabel='Per-gene β (log10 years per SD)', ylabel='−log10 P (Wald)', title='Exploratory · per gene')
    ax.legend(loc='upper left', fontsize=7.5)
    save(fig, 'fig6_per_gene')


def fig_template():
    path = RESULTS / 'template_sensitivity_species_M.csv'
    if not path.exists():
        return
    st = pd.read_csv(RESULTS / 'species_metrics.csv').set_index('binomial')
    m = pd.read_csv(path, index_col=0).join(st[['group']])
    info = json.loads((RESULTS / 'template_sensitivity.json').read_text())
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.7))
    m = m.reset_index()
    by_group(axes[0], m.zebrafish.to_numpy(), m.latimeria_chalumnae.to_numpy(), m)
    axes[0].set(xlabel='ΔΔG/difference · zebrafish template', ylabel='ΔΔG/difference · coelacanth template',
                title=f"All species ({info['genes_used']} genes)")
    axes[0].text(0.98, 0.03, f"r = {info['r_species_zebrafish_coelacanth']:.2f}", transform=axes[0].transAxes, ha='right', fontsize=7.5, color=INK2)
    rep = m[m.group != 'Amphibia']
    by_group(axes[1], rep.zebrafish.to_numpy(), rep.homo_sapiens.to_numpy(), rep)
    axes[1].set(xlabel='ΔΔG/difference · zebrafish template', ylabel='ΔΔG/difference · human template',
                title='Reptiles only')
    axes[1].text(0.98, 0.03, f"r = {info['r_species_reptiles_zebrafish_human']:.2f}", transform=axes[1].transAxes, ha='right', fontsize=7.5, color=INK2)
    axes[0].legend(loc='upper left', fontsize=7)
    fig.tight_layout()
    save(fig, 'figS1_templates')


def main():
    st = pd.read_csv(RESULTS / 'species_metrics.csv')
    s = json.loads((RESULTS / 'stats_summary.json').read_text())
    c = json.loads((RESULTS / 'predictor_controls_summary.json').read_text())
    fig_species(st)
    fig_controls(c)
    fig_primary(st, s)
    fig_competitive(s)
    fig_sensitivity(s)
    fig_genes()
    fig_template()
    print('figures written to', OUT)


if __name__ == '__main__':
    main()
