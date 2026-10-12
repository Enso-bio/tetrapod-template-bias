"""LaTeX tables and number macros for the manuscript, generated from the result files (no hand-typed values)."""
import json

import pandas as pd

from common import DATA, RESULTS, ROOT

OUT = ROOT / 'paper' / 'generated'
OUT.mkdir(parents=True, exist_ok=True)


def write_macro(name, rows):
    r"""Table bodies as macros: \input inside tabularx breaks booktabs rules."""
    body = '\n'.join(rows)
    (OUT / f'{name}.tex').write_text(f'\\newcommand{{\\{name}}}{{%\n{body}\n}}\n', encoding='utf-8')


def num(x, d=3):
    s = f'{x:.{d}f}'
    return s.replace('-', '$-$') if s.startswith('-') else s


def pval(p):
    return '$<$0.001' if p < 0.001 else (f'{p:.3f}' if p < 0.1 else f'{p:.2f}')


def macros():
    s = json.loads((RESULTS / 'stats_summary.json').read_text())
    o = json.loads((RESULTS / 'posthoc_monodelphis_domestica_template.json').read_text())
    c = json.loads((RESULTS / 'predictor_controls_summary.json').read_text())
    t = json.loads((RESULTS / 'template_sensitivity.json').read_text())
    sens = pd.read_csv(RESULTS / 'sensitivity_P1.csv').set_index('analysis')
    post = sens[sens.index.str.startswith('posterior_tree')]
    m = {
        'NSpecies': s['n_species'], 'NGenes': f"{s['n_genes']:,}", 'NGenAge': s['n_genage'],
        'NBackground': f"{s['n_background']:,}", 'NGenesOpossum': f"{o['n_genes']:,}",
        'POneBeta': num(s['P1']['beta']), 'POneLo': num(s['P1']['ci95_low']), 'POneHi': num(s['P1']['ci95_high']),
        'POneP': f"{s['P1']['p_boot']:.3f}", 'POneLambda': f"{s['P1']['lambda']:.2f}",
        'PTwoBeta': num(s['P2_selfcontained']['beta']), 'PTwoLo': num(s['P2_selfcontained']['ci95_low']),
        'PTwoHi': num(s['P2_selfcontained']['ci95_high']), 'PTwoP': f"{s['P2_selfcontained']['p_boot']:.3f}",
        'PTwoCompP': f"{s['P2_competitive']['p_empirical_two_sided']:.3f}",
        'PTwoNullMean': num(s['P2_competitive']['null_mean']), 'PTwoNullSD': num(s['P2_competitive']['null_sd']),
        'OpBeta': num(o['P1_posthoc']['beta']), 'OpLo': num(o['P1_posthoc']['ci95_low']), 'OpHi': num(o['P1_posthoc']['ci95_high']),
        'OpP': f"{o['P1_posthoc']['p_boot']:.3f}", 'OpCompP': f"{o['P2_competitive_posthoc']['p_empirical_two_sided']:.3f}",
        'OpDetect': f"{100 * (10 ** ((1.96 + 0.8416) * o['P1_posthoc']['se']) - 1):.0f}",
        'MainDetect': f"{s['P1_detectable_lifespan_percent_per_SD']:.0f}",
        'RBeta': num(s['secondary']['R_all']['beta']), 'RLo': num(s['secondary']['R_all']['ci95_low']),
        'RHi': num(s['secondary']['R_all']['ci95_high']), 'RP': f"{s['secondary']['R_all']['p_boot']:.3f}",
        'NP': f"{s['secondary']['N_all']['p_boot']:.2f}",
        'NNBeta': num(s['secondary']['nearly_neutral_M_on_log10mass']['beta']),
        'NNP': f"{s['secondary']['nearly_neutral_M_on_log10mass']['p_wald']:.2f}",
        'NNOpBeta': num(o['nearly_neutral_M_on_log10mass_posthoc']['beta']),
        'NNOpP': f"{o['nearly_neutral_M_on_log10mass_posthoc']['p_wald']:.2f}",
        'ExPrimBeta': num(sens.loc['exclude_primates', 'beta']), 'ExPrimLo': num(sens.loc['exclude_primates', 'ci95_low']),
        'ExPrimHi': num(sens.loc['exclude_primates', 'ci95_high']),
        'ExHumBeta': num(sens.loc['exclude_human', 'beta']), 'ExHumP': f"{sens.loc['exclude_human', 'p_wald']:.2f}",
        'AdjRBeta': num(sens.loc['adjust_log10_R', 'beta']), 'AdjRP': f"{sens.loc['adjust_log10_R', 'p_wald']:.2f}",
        'OLSP': pval(sens.loc['OLS_HC3', 'p_wald']),
        'PostBetaMin': num(post.beta.min()), 'PostBetaMax': num(post.beta.max()),
        'PostPMin': f"{post.p_wald.min():.2f}", 'PostPMax': f"{post.p_wald.max():.2f}",
        'MouseSiteR': f"{t['site_ddG_correlation_human_vs_mouse_template']:.2f}",
        'MouseSpeciesR': f"{t['species_M_correlation']:.2f}",
        'MouseBetaHuman': num(t['P1_like_human_template']['beta']), 'MouseBetaMouse': num(t['P1_like_mouse_template']['beta']),
        'COneObs': f"{c['C1_mean_observed_ddG']:.2f}", 'COneAlt': f"{c['C1_mean_alternative_ddG']:.2f}",
        'COneFrac': f"{100 * c['C1_fraction_genes_observed_below_alternatives']:.1f}",
        'CTwoBur': f"{c['C2_mean_abs_alt_buried']:.2f}", 'CTwoExp': f"{c['C2_mean_abs_alt_exposed']:.2f}",
        'CTwoFrac': f"{100 * c['C2_fraction_genes_buried_gt_exposed']:.1f}",
        'GenesTested': f"{s['per_gene']['tested']:,}", 'GenesPos': f"{100 * s['per_gene']['fraction_beta_positive']:.1f}",
    }
    e5_path = RESULTS / 'E5_confirmatory.json'
    if e5_path.exists():
        e = json.loads(e5_path.read_text())
        d = e['descriptive']
        m.update({
            'EGenes': f"{e['n_common_genes']:,}", 'EOffHum': num(e['primate_offset_SD']['human'], 2),
            'EOffOpo': num(e['primate_offset_SD']['opossum'], 2), 'EOffPla': num(e['primate_offset_SD']['platypus'], 2),
            'EaEst': num(e['H5a']['estimate'], 2), 'EaLo': num(e['H5a']['ci95'][0], 2), 'EaHi': num(e['H5a']['ci95'][1], 2),
            'EaVerdict': 'confirmed' if e['H5a']['confirmed'] else 'not confirmed',
            'EbEst': num(e['H5b']['estimate'], 2), 'EbLo': num(e['H5b']['ci95'][0], 2), 'EbHi': num(e['H5b']['ci95'][1], 2),
            'ErOpPl': f"{e['H5b']['r_opossum_platypus']:.2f}", 'ErHuPl': f"{e['H5b']['r_human_platypus']:.2f}",
            'ErHuOp': f"{e['H5b']['r_human_opossum']:.2f}",
            'EbVerdict': 'confirmed' if e['H5b']['confirmed'] else 'not confirmed',
            'EcBeta': num(e['H5c']['beta']), 'EcLo': num(e['H5c']['ci95_low']), 'EcHi': num(e['H5c']['ci95_high']),
            'EcP': f"{e['H5c']['p_boot']:.3f}", 'EcVerdict': 'replicated' if e['H5c']['replicated_null'] else 'not replicated',
            'EBetaHum': num(d['P1_slopes_common_sites_wald']['human']['beta']),
            'EBetaOpo': num(d['P1_slopes_common_sites_wald']['opossum']['beta']),
            'ERankHum': d['human_rank_of_57']['human'], 'ERankOpo': d['human_rank_of_57']['opossum'],
            'ERankPla': d['human_rank_of_57']['platypus'],
            'ECompP': f"{d['genage_competitive_platypus']['p_empirical_two_sided']:.3f}",
        })
    lines = ['% Generated by src/11_paper_tables.py from results/*.json. Do not edit by hand.']
    lines += [f'\\newcommand{{\\{k}}}{{{v}}}' for k, v in m.items()]
    (OUT / 'numbers.tex').write_text('\n'.join(lines) + '\n', encoding='utf-8')


def species_table():
    sp = pd.read_csv(DATA / 'species_final.csv').set_index('ensembl_name')
    st = pd.read_csv(RESULTS / 'species_metrics.csv')
    o = pd.read_csv(RESULTS / 'species_metrics_monodelphis_domestica.csv').set_index('ensembl_name')
    st['binomial'] = st.ensembl_name.map(sp.anage_binomial).str.replace('_', ' ').str.capitalize()
    st['M_op'] = st.ensembl_name.map(o.M_all_outgroup)
    groups = {'corta': 'S', 'intermedia': 'I', 'larga': 'L'}
    rows = []
    for r in st.sort_values(['order', 'binomial']).itertuples():
        rows.append(f"\\emph{{{r.binomial}}} & {r.order} & {r.max_lifespan_yr:g} & {r.adult_mass_g:,.0f} & "
                    f"{r.LQ:.2f} & {groups[r.LQ_group]} & {r.M_all:.3f} & {r.M_op:.3f} & {r.R_all:.4f} \\\\")
    write_macro('SpeciesRows', rows)


def primary_table():
    s = json.loads((RESULTS / 'stats_summary.json').read_text())
    o = json.loads((RESULTS / 'posthoc_monodelphis_domestica_template.json').read_text())
    def row(name, r, p, extra=''):
        return (f"{name} & {num(r['beta'])} & [{num(r['ci95_low'])}, {num(r['ci95_high'])}] & "
                f"{r['lambda']:.2f} & {p}{extra} \\\\")
    rows = [
        row('P1: all genes, human template', s['P1'], f"{s['P1']['p_boot']:.3f}"),
        row('P2: GenAge genes, human template', s['P2_selfcontained'], f"{s['P2_selfcontained']['p_boot']:.3f}"),
        f"P2 competitive: GenAge vs 1,000 background sets & -- & -- & -- & {s['P2_competitive']['p_empirical_two_sided']:.3f}$^{{a}}$ \\\\",
        '\\midrule',
        row('Divergence from ancestor per site (secondary)', s['secondary']['R_all'], f"{s['secondary']['R_all']['p_boot']:.3f}"),
        row('Net $\\Delta\\Delta G$ per site (secondary)', s['secondary']['N_all'], f"{s['secondary']['N_all']['p_boot']:.3f}"),
        '\\midrule',
        row('P1, opossum template (post hoc)', o['P1_posthoc'], f"{o['P1_posthoc']['p_boot']:.3f}"),
        row('P2, opossum template (post hoc)', o['P2_selfcontained_posthoc'], f"{o['P2_selfcontained_posthoc']['p_boot']:.3f}"),
        f"P2 competitive, opossum template (post hoc) & -- & -- & -- & {o['P2_competitive_posthoc']['p_empirical_two_sided']:.3f}$^{{a}}$ \\\\",
    ]
    e5_path = RESULTS / 'E5_confirmatory.json'
    if e5_path.exists():
        e = json.loads(e5_path.read_text())
        c = e['descriptive']['genage_competitive_platypus']
        rows += ['\\midrule',
                 row('P1, platypus template, common sites (confirmatory, H5c)', e['H5c'], f"{e['H5c']['p_boot']:.3f}"),
                 f"P2 competitive, platypus template, common sites & -- & -- & -- & {c['p_empirical_two_sided']:.3f}$^{{a}}$ \\\\"]
    write_macro('PrimaryRows', rows)


def sensitivity_table():
    sens = pd.read_csv(RESULTS / 'sensitivity_P1.csv')
    names = {'lambda_fixed_1': '$\\lambda$ fixed at 1', 'variant_no_mask': 'No window mask',
             'variant_plddt90': 'Sites with pLDDT $\\geq$ 90', 'variant_buried': 'Buried sites only',
             'variant_exposed': 'Exposed sites only', 'exclude_primates': 'Primates excluded',
             'exclude_human': '\\emph{Homo sapiens} excluded', 'exclude_bats': 'Bats excluded',
             'high_quality_only': 'AnAge quality ``high\'\' only', 'exclude_rodents': 'Rodents excluded',
             'adjust_log10_R': 'Adjusted for log$_{10}$ divergence from ancestor', 'background_genes_only': 'Background genes only',
             'response_log10_LQ_no_mass': 'Response log$_{10}$ LQ (no mass term)', 'OLS_HC3': 'OLS with HC3 errors (no phylogeny)'}
    rows = []
    for _, r in sens[~sens.analysis.str.startswith('posterior')].iterrows():
        lam = '--' if pd.isna(r['lambda']) else f"{r['lambda']:.2f}"
        rows.append(f"{names.get(r.analysis, r.analysis)} & {int(r.n)} & {num(r.beta)} & "
                    f"[{num(r.ci95_low)}, {num(r.ci95_high)}] & {lam} & {pval(r.p_wald)} \\\\")
    post = sens[sens.analysis.str.startswith('posterior')]
    rows.append(f"20 posterior trees (range) & 57 & {num(post.beta.min())} to {num(post.beta.max())} & -- & -- & "
                f"{post.p_wald.min():.2f}--{post.p_wald.max():.2f} \\\\")
    write_macro('SensitivityRows', rows)


if __name__ == '__main__':
    macros()
    species_table()
    primary_table()
    sensitivity_table()
    print('written to', OUT)
