"""Publication figures (PNG + PDF + SVG) in Spanish (figures/) and English (paper/figures/).

Usage: python 09_figures.py [es|en|all]   (default: all)
"""
import json
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import FIGURES, RESULTS, ROOT

INK, INK2, MUTED, GRID, AXIS = '#0b0b0b', '#52514e', '#898781', '#e1e0d9', '#c3c2b7'
BLUE, ORANGE = '#2a78d6', '#eb6834'
ORDINAL = {'corta': '#86b6ef', 'intermedia': '#2a78d6', 'larga': '#104281'}
plt.rcParams.update({'font.family': 'sans-serif', 'font.sans-serif': ['Segoe UI', 'DejaVu Sans'], 'font.size': 9,
                     'axes.edgecolor': AXIS, 'axes.labelcolor': INK2, 'xtick.color': MUTED, 'ytick.color': MUTED,
                     'axes.titlesize': 10, 'axes.titleweight': 'bold', 'axes.titlecolor': INK,
                     'axes.spines.top': False, 'axes.spines.right': False, 'axes.grid': True, 'grid.color': GRID,
                     'grid.linewidth': 0.6, 'legend.frameon': False, 'figure.dpi': 150, 'savefig.dpi': 300,
                     'savefig.bbox': 'tight', 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
LABELS = {'homo_sapiens': 'Homo sapiens', 'heterocephalus_glaber_female': 'H. glaber', 'mus_musculus': 'Mus musculus',
          'myotis_lucifugus': 'M. lucifugus', 'balaenoptera_musculus': 'B. musculus', 'physeter_catodon': 'P. catodon',
          'loxodonta_africana': 'L. africana', 'rattus_norvegicus': 'R. norvegicus'}

TEXT = {
    'es': {
        'out': FIGURES,
        'group': {'corta': 'corta', 'intermedia': 'intermedia', 'larga': 'larga'},
        'group_label': 'Longevidad relativa {g} (n={n})', 'expected': 'Esperado por masa (LQ = 1)',
        'mass': 'Masa adulta (g)', 'lifespan': 'Vida máxima registrada (años)', 'species_title': '{n} mamíferos placentarios analizados',
        'obs': 'Sustituciones\nobservadas', 'alt': 'Alternativas\nposibles', 'ddg_mean': 'ΔΔG medio predicho (kcal/mol)',
        'c1': 'C1 · selección purificadora', 'buried': 'Sitios\nenterrados', 'exposed': 'Sitios\nexpuestos',
        'abs_alt': '|ΔΔG| medio de alternativas', 'c2': 'C2 · enterramiento',
        'c3rows': ['Hélice → Pro', 'Hélice → Gly', 'Hélice (media)', 'Lámina → Pro', 'Lámina → Gly', 'Lámina (media)'],
        'ddg_mean_short': 'ΔΔG medio predicho', 'c3': 'C3 · Pro/Gly en estructura',
        'others': 'Otros placentarios (n={n})', 'primates': 'Primates (n={n})',
        'x_main': 'ΔΔG medio por sustitución (DE, todos los genes; plantilla humana)',
        'y_lq': 'log10 LQ (longevidad relativa a la masa)', 'main_title': 'P1 · estabilidad predicha y longevidad',
        'main_txt': 'PGLS: β = {b:.3f} log10 años/DE\nIC95 {lo:.3f} a {hi:.3f}\np bootstrap = {p:.3f} · λ = {lam:.2f} · n = {n}',
        'human_tpl': 'Plantilla humana (P1)', 'opossum_tpl': 'Plantilla de zarigüeya (post hoc)',
        'x_sd': 'ΔΔG medio por sustitución (DE)', 'beta_txt': 'β = {b:.3f} (IC95 {lo:.3f} a {hi:.3f})\np bootstrap = {p:.3f}',
        'genage_txt': 'GenAge (n={n})\np empírico = {p:.3f}', 'x_comp': 'β de ΔΔG por sustitución (log10 años/DE)',
        'y_comp': 'Conjuntos aleatorios de fondo', 'comp_title': 'P2 · GenAge frente a {n} conjuntos de fondo',
        'sens': {'lambda_fixed_1': 'λ = 1 fijo', 'variant_no_mask': 'Sin filtro de ventanas', 'variant_plddt90': 'Solo pLDDT ≥ 90',
                 'variant_buried': 'Solo sitios enterrados', 'variant_exposed': 'Solo sitios expuestos',
                 'exclude_primates': 'Sin primates', 'exclude_human': 'Sin humano', 'exclude_bats': 'Sin murciélagos',
                 'exclude_rodents': 'Sin roedores', 'high_quality_only': 'Solo AnAge "high"',
                 'adjust_log10_R': 'Ajuste por tasa de sustitución', 'background_genes_only': 'Solo genes de fondo',
                 'response_log10_LQ_no_mass': 'Respuesta log10 LQ', 'OLS_HC3': 'MCO, errores HC3'},
        'primary': 'Principal (P1)', 'posterior': '20 árboles posteriores (rango β)',
        'x_sens': 'β (log10 años por DE) e IC95 aproximado', 'sens_title': 'Sensibilidad de P1',
        'bg': 'Fondo (n={n})', 'ga': 'GenAge y familias previas (n={n})', 'old': 'Familias del estudio anterior (n={n})',
        'fdr': ' FDR 5 %', 'x_gene': 'β por gen (log10 años por DE)', 'gene_title': 'Exploratorio · gen a gen',
        'x_tpl': 'ΔΔG/sustitución · plantilla humana', 'y_tpl': 'ΔΔG/sustitución · plantilla de ratón',
        'tpl_title': 'Plantilla estructural ({n} genes)', 'tpl_txt': 'r especies = {a:.3f}\nr sitios = {b:.3f}',
        'e5_templates': ['Humana', 'Zarigüeya', 'Ornitorrinco'], 'e5_offset': 'Desplazamiento de primates (DE) e IC95',
        'e5_title_a': 'E5 · sesgo hacia primates', 'e5_x': 'ΔΔG/sustitución · zarigüeya',
        'e5_y': 'ΔΔG/sustitución · ornitorrinco', 'e5_title_b': 'E5 · dos plantillas externas',
        'e5_r': 'r zarigüeya–ornitorrinco = {a:.2f}\nr humana–ornitorrinco = {b:.2f}',
        'files': {'confirm': 'fig8_confirmatoria_E5', 'species': 'fig1_especies', 'controls': 'fig2_controles_predictor', 'main': 'fig3_resultado_principal',
                  'comp': 'fig4_competitiva_genage', 'sens': 'fig5_sensibilidad', 'genes': 'fig6_por_gen',
                  'outgroup': 'fig7_sesgo_plantilla', 'template': 'figS1_plantilla'}},
    'en': {
        'out': ROOT / 'paper' / 'figures',
        'group': {'corta': 'short', 'intermedia': 'intermediate', 'larga': 'long'},
        'group_label': '{g} relative longevity (n={n})', 'expected': 'Expected from mass (LQ = 1)',
        'mass': 'Adult mass (g)', 'lifespan': 'Maximum recorded lifespan (years)', 'species_title': '{n} placental mammals analysed',
        'obs': 'Observed\ndifferences', 'alt': 'Possible\nalternatives', 'ddg_mean': 'Mean predicted ΔΔG (kcal/mol)',
        'c1': 'C1 · purifying selection', 'buried': 'Buried\nsites', 'exposed': 'Exposed\nsites',
        'abs_alt': 'Mean |ΔΔG| of alternatives', 'c2': 'C2 · burial',
        'c3rows': ['Helix → Pro', 'Helix → Gly', 'Helix (mean)', 'Strand → Pro', 'Strand → Gly', 'Strand (mean)'],
        'ddg_mean_short': 'Mean predicted ΔΔG', 'c3': 'C3 · Pro/Gly in structure',
        'others': 'Other placentals (n={n})', 'primates': 'Primates (n={n})',
        'x_main': 'Mean ΔΔG per difference (SD; all genes; human template)',
        'y_lq': 'log10 LQ (lifespan relative to mass)', 'main_title': 'P1 · predicted stability and lifespan',
        'main_txt': 'PGLS: β = {b:.3f} log10 years/SD\n95% CI {lo:.3f} to {hi:.3f}\nbootstrap P = {p:.3f} · λ = {lam:.2f} · n = {n}',
        'human_tpl': 'Human template (P1)', 'opossum_tpl': 'Opossum template (post hoc)',
        'x_sd': 'Mean ΔΔG per difference (SD)', 'beta_txt': 'β = {b:.3f} (95% CI {lo:.3f} to {hi:.3f})\nbootstrap P = {p:.3f}',
        'genage_txt': 'GenAge (n={n})\nempirical P = {p:.3f}', 'x_comp': 'β of ΔΔG per difference (log10 years/SD)',
        'y_comp': 'Random background sets', 'comp_title': 'P2 · GenAge vs {n} background sets',
        'sens': {'lambda_fixed_1': 'λ fixed at 1', 'variant_no_mask': 'No window mask', 'variant_plddt90': 'pLDDT ≥ 90 only',
                 'variant_buried': 'Buried sites only', 'variant_exposed': 'Exposed sites only',
                 'exclude_primates': 'Without primates', 'exclude_human': 'Without human', 'exclude_bats': 'Without bats',
                 'exclude_rodents': 'Without rodents', 'high_quality_only': 'AnAge "high" only',
                 'adjust_log10_R': 'Adjusted for divergence from ancestor', 'background_genes_only': 'Background genes only',
                 'response_log10_LQ_no_mass': 'Response log10 LQ', 'OLS_HC3': 'OLS, HC3 errors'},
        'primary': 'Primary (P1)', 'posterior': '20 posterior trees (β range)',
        'x_sens': 'β (log10 years per SD) and approximate 95% CI', 'sens_title': 'Sensitivity of P1',
        'bg': 'Background (n={n})', 'ga': 'GenAge and previous panel (n={n})', 'old': 'Previous-study genes (n={n})',
        'fdr': ' FDR 5%', 'x_gene': 'Per-gene β (log10 years per SD)', 'gene_title': 'Exploratory · per gene',
        'x_tpl': 'ΔΔG/difference · human template', 'y_tpl': 'ΔΔG/difference · mouse template',
        'tpl_title': 'Structural template ({n} genes)', 'tpl_txt': 'r species = {a:.3f}\nr sites = {b:.3f}',
        'e5_templates': ['Human', 'Opossum', 'Platypus'], 'e5_offset': 'Primate offset (SD) and 95% CI',
        'e5_title_a': 'E5 · bias towards primates', 'e5_x': 'ΔΔG/difference · opossum template',
        'e5_y': 'ΔΔG/difference · platypus template', 'e5_title_b': 'E5 · two outgroup templates',
        'e5_r': 'r opossum–platypus = {a:.2f}\nr human–platypus = {b:.2f}',
        'files': {'confirm': 'fig8_confirmatory_E5', 'species': 'fig1_species', 'controls': 'fig2_controls', 'main': 'fig3_primary',
                  'comp': 'fig4_competitive', 'sens': 'fig5_sensitivity', 'genes': 'fig6_per_gene',
                  'outgroup': 'fig7_template_bias', 'template': 'figS1_mouse_template'}},
}
s_main = json.loads((RESULTS / 'stats_summary.json').read_text()) if (RESULTS / 'stats_summary.json').exists() else {}


class Figures:
    def __init__(self, lang):
        self.t = TEXT[lang]
        self.out = self.t['out']
        self.out.mkdir(parents=True, exist_ok=True)

    def save(self, fig, key):
        for ext in ('png', 'pdf', 'svg'):
            fig.savefig(self.out / f"{self.t['files'][key]}.{ext}")
        plt.close(fig)

    def species(self, st):
        t = self.t
        fig, ax = plt.subplots(figsize=(6.2, 4.2))
        for g in ['corta', 'intermedia', 'larga']:
            d = st[st.LQ_group == g]
            lab = t['group_label'].format(g=t['group'][g], n=len(d))
            ax.scatter(d.adult_mass_g, d.max_lifespan_yr, s=34, color=ORDINAL[g], edgecolor='white', linewidth=1.2,
                       label=lab[0].upper() + lab[1:], zorder=3)
        offsets = {'heterocephalus_glaber_female': (5, -10)}
        for k, lab in LABELS.items():
            r = st[st.ensembl_name == k]
            if len(r):
                ax.annotate(lab, (r.adult_mass_g.iloc[0], r.max_lifespan_yr.iloc[0]), xytext=offsets.get(k, (5, 3)),
                            textcoords='offset points', fontsize=7.5, color=INK2, style='italic')
        m = np.logspace(np.log10(st.adult_mass_g.min() / 2), np.log10(st.adult_mass_g.max() * 2), 50)
        ax.plot(m, 4.88 * m ** 0.153, color=MUTED, lw=1, ls='--', zorder=2, label=t['expected'])
        ax.set(xscale='log', yscale='log', xlabel=t['mass'], ylabel=t['lifespan'], title=t['species_title'].format(n=len(st)))
        ax.legend(loc='upper left', fontsize=7.5)
        self.save(fig, 'species')

    def controls(self, c):
        t = self.t
        fig, axes = plt.subplots(1, 3, figsize=(8.4, 2.9))
        axes[0].barh([1, 0], [c['C1_mean_observed_ddG'], c['C1_mean_alternative_ddG']], color=[BLUE, AXIS], height=0.5)
        axes[0].set_yticks([1, 0], [t['obs'], t['alt']])
        axes[0].set(xlabel=t['ddg_mean'], title=t['c1'])
        axes[1].barh([1, 0], [c['C2_mean_abs_alt_buried'], c['C2_mean_abs_alt_exposed']], color=[BLUE, AXIS], height=0.5)
        axes[1].set_yticks([1, 0], [t['buried'], t['exposed']])
        axes[1].set(xlabel=t['abs_alt'], title=t['c2'])
        vals = [c['C3_helix_region_to_P'], c['C3_helix_region_to_G'], c['C3_helix_region_mean_alt'],
                c['C3_strand_region_to_P'], c['C3_strand_region_to_G'], c['C3_strand_region_mean_alt']]
        axes[2].barh(range(6)[::-1], vals, height=0.6, color=[BLUE, BLUE, AXIS, BLUE, BLUE, AXIS])
        axes[2].set_yticks(range(6)[::-1], t['c3rows'])
        axes[2].set(xlabel=t['ddg_mean_short'], title=t['c3'])
        for a in axes:
            a.grid(axis='y', visible=False)
        fig.tight_layout()
        self.save(fig, 'controls')

    def scatter_by_primate(self, ax, x, y, st):
        prim = (st.order == 'Primates').to_numpy()
        ax.scatter(x[~prim], y[~prim], s=30, color=BLUE, edgecolor='white', linewidth=1.2, zorder=3,
                   label=self.t['others'].format(n=(~prim).sum()))
        ax.scatter(x[prim], y[prim], s=30, color=ORANGE, edgecolor='white', linewidth=1.2, zorder=3,
                   label=self.t['primates'].format(n=prim.sum()))
        for k, lab, off in [('homo_sapiens', 'Homo sapiens', (6, -10)), ('mus_musculus', 'Mus musculus', (6, 2)),
                            ('heterocephalus_glaber_female', 'H. glaber', (6, -9)), ('myotis_lucifugus', 'M. lucifugus', (6, 2))]:
            i = np.where(st.ensembl_name.to_numpy() == k)[0]
            if len(i):
                ax.annotate(lab, (x[i[0]], y[i[0]]), xytext=off, textcoords='offset points', fontsize=7.5,
                            color=INK2, style='italic')

    def main_fig(self, st, s):
        t, p1 = self.t, s['P1']
        x = ((st.M_all - st.M_all.mean()) / st.M_all.std(ddof=1)).to_numpy()
        fig, ax = plt.subplots(figsize=(5.4, 4.0))
        self.scatter_by_primate(ax, x, np.log10(st.LQ.to_numpy()), st)
        ax.set(xlabel=t['x_main'], ylabel=t['y_lq'], title=t['main_title'])
        ax.text(0.02, 0.03, t['main_txt'].format(b=p1['beta'], lo=p1['ci95_low'], hi=p1['ci95_high'], p=p1['p_boot'],
                                                 lam=p1['lambda'], n=p1['n']),
                transform=ax.transAxes, va='bottom', fontsize=7.5, color=INK2)
        ax.legend(loc='upper center', fontsize=7.5)
        self.save(fig, 'main')

    def outgroup(self):
        path = RESULTS / 'species_metrics_monodelphis_domestica.csv'
        if not path.exists():
            return
        t = self.t
        o = pd.read_csv(path)
        info = json.loads((RESULTS / 'posthoc_monodelphis_domestica_template.json').read_text())
        zf = lambda v: ((v - v.mean()) / v.std(ddof=1)).to_numpy()
        fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.8), sharey=True)
        for ax, col, title, r, (tx, ha) in [(axes[0], 'M_all_human_template', t['human_tpl'], s_main['P1'], (0.02, 'left')),
                                            (axes[1], 'M_all_outgroup', t['opossum_tpl'], info['P1_posthoc'], (0.98, 'right'))]:
            self.scatter_by_primate(ax, zf(o[col]), np.log10(o.LQ.to_numpy()), o)
            ax.set(xlabel=t['x_sd'], title=title)
            ax.text(tx, 0.03, t['beta_txt'].format(b=r['beta'], lo=r['ci95_low'], hi=r['ci95_high'], p=r['p_boot']),
                    transform=ax.transAxes, va='bottom', ha=ha, fontsize=7.5, color=INK2)
        axes[0].set_ylabel('log10 LQ')
        axes[1].legend(loc='upper right', fontsize=7.5)
        fig.tight_layout()
        self.save(fig, 'outgroup')

    def competitive(self, s):
        t, c = self.t, s['P2_competitive']
        null = np.load(RESULTS / 'P2_null_betas.npy')
        fig, ax = plt.subplots(figsize=(5.4, 3.2))
        ax.hist(null, bins=40, color=AXIS, edgecolor='white', linewidth=0.8)
        ax.axvline(c['beta_genage'], color=ORANGE, lw=2)
        ax.annotate(t['genage_txt'].format(n=c['set_size'], p=c['p_empirical_two_sided']),
                    (c['beta_genage'], ax.get_ylim()[1] * 0.9), xytext=(6, 0), textcoords='offset points',
                    fontsize=7.5, color=INK2)
        ax.set(xlabel=t['x_comp'], ylabel=t['y_comp'], title=t['comp_title'].format(n=c['n_null_sets']))
        ax.grid(axis='x', visible=False)
        self.save(fig, 'comp')

    def sensitivity(self, s):
        t = self.t
        sens = pd.read_csv(RESULTS / 'sensitivity_P1.csv')
        post = sens[sens.analysis.str.startswith('posterior_tree')]
        other = sens[~sens.analysis.str.startswith('posterior_tree') & sens.beta.notna()]
        rows = [(t['primary'], s['P1']['beta'], s['P1']['ci95_low'], s['P1']['ci95_high'])]
        rows += [(t['sens'].get(r.analysis, r.analysis), r.beta, r.ci95_low, r.ci95_high) for r in other.itertuples()]
        rows.append((t['posterior'], post.beta.median(), post.beta.min(), post.beta.max()))
        fig, ax = plt.subplots(figsize=(5.6, 0.28 * len(rows) + 0.8))
        ys = np.arange(len(rows))[::-1]
        for k, (y, (n, b, lo, hi)) in enumerate(zip(ys, rows)):
            col = ORANGE if k == 0 else BLUE
            ax.plot([lo, hi], [y, y], color=col, lw=2, solid_capstyle='round')
            ax.scatter([b], [y], color=col, s=26, zorder=3, edgecolor='white', linewidth=1)
        ax.axvline(0, color=MUTED, lw=1, ls='--')
        ax.set_yticks(ys, [r[0] for r in rows])
        ax.set(xlabel=t['x_sens'], title=t['sens_title'])
        ax.grid(axis='y', visible=False)
        self.save(fig, 'sens')

    def genes(self):
        t = self.t
        g = pd.read_csv(RESULTS / 'per_gene_exploratory.csv')
        fig, ax = plt.subplots(figsize=(5.6, 3.8))
        bg, ga, old = g[g.gene_set == 'background'], g[g.gene_set != 'background'], g[g.old_family]
        ax.scatter(bg.beta, -np.log10(bg.p_wald), s=9, color=AXIS, label=t['bg'].format(n=len(bg)), zorder=2)
        ax.scatter(ga.beta, -np.log10(ga.p_wald), s=12, color=ORANGE, label=t['ga'].format(n=len(ga)), zorder=3)
        ax.scatter(old.beta, -np.log10(old.p_wald), s=40, facecolor='none', edgecolor=INK, linewidth=1.1, zorder=4,
                   label=t['old'].format(n=len(old)))
        if (g.q_bh <= 0.05).any():
            thr = g.loc[g.q_bh <= 0.05, 'p_wald'].max()
            ax.axhline(-np.log10(thr), color=MUTED, lw=1, ls='--')
            ax.text(ax.get_xlim()[0], -np.log10(thr), t['fdr'], va='bottom', fontsize=7, color=MUTED)
        ax.set(xlabel=t['x_gene'], ylabel='−log10 P (Wald)', title=t['gene_title'])
        ax.legend(loc='upper left', bbox_to_anchor=(0, 0.9), fontsize=7.5)
        self.save(fig, 'genes')

    def confirm(self):
        path = RESULTS / 'E5_confirmatory.json'
        if not path.exists():
            return
        t = self.t
        e = json.loads(path.read_text())
        m = pd.read_csv(RESULTS / 'E5_species_metrics_common_sites.csv')
        fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.5), gridspec_kw={'width_ratios': [1, 1.25]})
        ax = axes[0]
        keys = [('human', 'off_h'), ('opossum', 'off_o'), ('platypus', 'off_p')]
        for k, (name, bk) in enumerate(keys):
            lo, hi = e['descriptive']['bootstrap_offset_ci95'][bk]
            col = ORANGE if name == 'human' else BLUE
            ax.plot([lo, hi], [2 - k, 2 - k], color=col, lw=2, solid_capstyle='round')
            ax.scatter([e['primate_offset_SD'][name]], [2 - k], color=col, s=30, zorder=3, edgecolor='white', linewidth=1)
        ax.axvline(0, color=MUTED, lw=1, ls='--')
        ax.set_yticks([2, 1, 0], t['e5_templates'])
        ax.set(xlabel=t['e5_offset'], title=t['e5_title_a'])
        ax.grid(axis='y', visible=False)
        ax = axes[1]
        prim = (m.order == 'Primates').to_numpy()
        ax.scatter(m.M_opossum[~prim], m.M_platypus[~prim], s=28, color=BLUE, edgecolor='white', linewidth=1.2, zorder=3,
                   label=t['others'].format(n=(~prim).sum()))
        ax.scatter(m.M_opossum[prim], m.M_platypus[prim], s=28, color=ORANGE, edgecolor='white', linewidth=1.2, zorder=3,
                   label=t['primates'].format(n=prim.sum()))
        h = m[m.ensembl_name == 'homo_sapiens']
        ax.annotate('Homo sapiens', (h.M_opossum.iloc[0], h.M_platypus.iloc[0]), xytext=(6, -10), textcoords='offset points',
                    fontsize=7.5, color=INK2, style='italic')
        ax.set(xlabel=t['e5_x'], ylabel=t['e5_y'], title=t['e5_title_b'])
        ax.text(0.98, 0.03, t['e5_r'].format(a=e['H5b']['r_opossum_platypus'], b=e['H5b']['r_human_platypus']),
                transform=ax.transAxes, ha='right', va='bottom', fontsize=7.5, color=INK2)
        ax.legend(loc='upper left', fontsize=7.5)
        fig.tight_layout()
        self.save(fig, 'confirm')

    def template(self):
        path = RESULTS / 'template_sensitivity_species_M.csv'
        if not path.exists():
            return
        t = self.t
        d = pd.read_csv(path)
        info = json.loads((RESULTS / 'template_sensitivity.json').read_text())
        fig, ax = plt.subplots(figsize=(4.2, 3.8))
        ax.scatter(d.human, d.mouse, s=26, color=BLUE, edgecolor='white', linewidth=1.2, zorder=3)
        lo, hi = min(d.human.min(), d.mouse.min()), max(d.human.max(), d.mouse.max())
        ax.plot([lo, hi], [lo, hi], color=MUTED, lw=1, ls='--')
        ax.set(xlabel=t['x_tpl'], ylabel=t['y_tpl'], title=t['tpl_title'].format(n=info['genes_used']))
        ax.text(0.02, 0.98, t['tpl_txt'].format(a=info['species_M_correlation'],
                                                b=info['site_ddG_correlation_human_vs_mouse_template']),
                transform=ax.transAxes, va='top', fontsize=7.5, color=INK2)
        self.save(fig, 'template')


def main():
    langs = ['es', 'en'] if len(sys.argv) < 2 or sys.argv[1] == 'all' else [sys.argv[1]]
    st = pd.read_csv(RESULTS / 'species_metrics.csv')
    c = json.loads((RESULTS / 'predictor_controls_summary.json').read_text())
    for lang in langs:
        f = Figures(lang)
        f.species(st)
        f.controls(c)
        f.main_fig(st, s_main)
        f.competitive(s_main)
        f.sensitivity(s_main)
        f.genes()
        f.template()
        f.outgroup()
        f.confirm()
        print('figures written to', f.out)


if __name__ == '__main__':
    main()
