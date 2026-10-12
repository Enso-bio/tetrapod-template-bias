"""LaTeX number macros and table bodies for the manuscript, generated from result files (no hand-typed values)."""
import json

import numpy as np
import pandas as pd

from common import DATA, RESULTS, ROOT

OUT = ROOT / 'paper' / 'generated'
OUT.mkdir(parents=True, exist_ok=True)
GROUP_NAMES = {'Testudines': 'Turtles', 'Crocodylia': 'Crocodilians', 'Lepidosauria': 'Squamates', 'Amphibia': 'Amphibians'}


def num(x, d=3):
    s = f'{x:.{d}f}'
    return s.replace('-', '$-$') if s.startswith('-') else s


def pval(p):
    return '$<$0.001' if p < 0.001 else (f'{p:.3f}' if p < 0.1 else f'{p:.2f}')


def write_macro(name, rows):
    body = '\n'.join(rows)
    (OUT / f'{name}.tex').write_text(f'\\newcommand{{\\{name}}}{{%\n{body}\n}}\n', encoding='utf-8')


def main():
    s = json.loads((RESULTS / 'stats_summary.json').read_text())
    c = json.loads((RESULTS / 'predictor_controls_summary.json').read_text())
    t = json.loads((RESULTS / 'template_sensitivity.json').read_text())
    sens = pd.read_csv(RESULTS / 'sensitivity_P1.csv').set_index('analysis')
    st = pd.read_csv(RESULTS / 'species_metrics.csv')
    met = pd.read_csv(RESULTS / 'species_gene_metrics.csv.gz')
    met = met[met.species.isin(st.binomial)]
    asr = pd.read_csv(RESULTS / 'asr_summary.csv')
    cand = pd.read_csv(DATA / 'species_candidates.csv')
    tt = json.loads((DATA / 'timetree_upgma_summary.json').read_text())
    sd = st.M_all.std(ddof=1)
    grp = lambda g: st.group == g
    m = {
        'NSpecies': s['n_species'], 'NTurtles': s['n_by_group'].get('Testudines', 0),
        'NCrocs': s['n_by_group'].get('Crocodylia', 0), 'NSquamates': s['n_by_group'].get('Lepidosauria', 0),
        'NAmphibians': s['n_by_group'].get('Amphibia', 0), 'NCandidates': int(cand.candidate.sum()),
        'NGenes': f"{s['n_genes']:,}", 'NGenAge': s['n_genage'], 'NBackground': f"{s['n_background']:,}",
        'NPairs': f"{len(met):,}", 'NSitesM': f"{met.L_main.sum() / 1e6:.1f}", 'NSubsM': f"{met.K_main.sum() / 1e6:.2f}",
        'MaskPct': f"{100 * met.masked_sites.sum() / (met.L_main.sum() + met.masked_sites.sum()):.1f}",
        'AsrPost': f"{asr.mean_max_posterior.mean():.3f}", 'RootAge': f"{tt['root_age_Ma']:.0f}",
        'POneBeta': num(s['P1']['beta']), 'POneLo': num(s['P1']['ci95_low']), 'POneHi': num(s['P1']['ci95_high']),
        'POneP': f"{s['P1']['p_boot']:.3f}", 'POneLambda': f"{s['P1']['lambda']:.2f}",
        'PTwoBeta': num(s['P2_selfcontained']['beta']), 'PTwoP': f"{s['P2_selfcontained']['p_boot']:.3f}",
        'PTwoCompP': f"{s['P2_competitive']['p_empirical_two_sided']:.3f}",
        'RBeta': num(s['secondary']['R']['beta']), 'RLo': num(s['secondary']['R']['ci95_low']),
        'RHi': num(s['secondary']['R']['ci95_high']), 'RP': f"{s['secondary']['R']['p_boot']:.3f}",
        'NNBeta': num(s['secondary']['nearly_neutral_M_on_log10mass']['beta']),
        'NNP': f"{s['secondary']['nearly_neutral_M_on_log10mass']['p_wald']:.2f}",
        'Detect': f"{s['P1_detectable_lifespan_percent_per_SD']:.0f}",
        'COneObs': f"{c['C1_mean_observed_ddG']:.2f}", 'COneAlt': f"{c['C1_mean_alternative_ddG']:.2f}",
        'COneFrac': f"{100 * c['C1_fraction_genes_observed_below_alternatives']:.1f}",
        'CTwoBur': f"{c['C2_mean_abs_alt_buried']:.2f}", 'CTwoExp': f"{c['C2_mean_abs_alt_exposed']:.2f}",
        'CTwoFrac': f"{100 * c['C2_fraction_genes_buried_gt_exposed']:.1f}",
        'GenesTested': f"{s['per_gene']['tested']:,}", 'GenesNominal': s['per_gene']['nominal_p_lt_0.05'],
        'GenesExpected': f"{0.05 * s['per_gene']['tested']:.0f}", 'GenesMinQ': f"{s['per_gene']['min_q']:.2f}",
        'GenesPos': f"{100 * s['per_gene']['fraction_beta_positive']:.1f}",
        'TempCov': f"{100 * sens.loc['coverage_body_temp_c_reptiles', 'coverage']:.0f}",
        'TplGenes': t['genes_used'], 'TplRCoel': f"{t['r_species_zebrafish_coelacanth']:.2f}",
        'TplBetaZf': num(t['P1_like_all_species_zebrafish']['beta']),
        'TplBetaCoel': num(t['P1_like_all_species_latimeria_chalumnae']['beta']),
        'TplRepZfBeta': num(t['P1_like_reptiles_only_zebrafish']['beta']),
        'TplRepZfP': f"{t['P1_like_reptiles_only_zebrafish']['p_wald']:.2f}",
        'TplRepHuBeta': num(t['P1_like_reptiles_only_homo_sapiens']['beta']),
        'TplRepHuP': f"{t['P1_like_reptiles_only_homo_sapiens']['p_wald']:.3f}",
        'TplRepHuLambda': f"{t['P1_like_reptiles_only_homo_sapiens']['lambda']:.2f}",
        'TplRRepHu': f"{t['r_species_reptiles_zebrafish_human']:.2f}",
        'AmphOffset': num((st.M_all[grp('Amphibia')].mean() - st.M_all[~grp('Amphibia')].mean()) / sd, 2),
        'TurtSqOffset': num((st.M_all[grp('Testudines')].mean() - st.M_all[grp('Lepidosauria')].mean()) / sd, 2),
    }
    for key in ['reptiles_only', 'amphibians_only', 'without_turtles', 'squamates_only', 'anage_mass_only',
                'no_mass_covariate', 'adjust_log10_R']:
        name = ''.join(w.capitalize() for w in key.replace('log10', 'log').replace('_', ' ').split())
        name = ''.join(ch for ch in name if ch.isalpha())
        m[f'S{name}Beta'] = num(sens.loc[key, 'beta'])
        m[f'S{name}P'] = f"{sens.loc[key, 'p_wald']:.2f}"
        m[f'S{name}N'] = int(sens.loc[key, 'n'])
    for g in GROUP_NAMES:
        m[f'Life{GROUP_NAMES[g]}'] = f"{st.loc[grp(g), 'max_lifespan_yr'].median():g}"
        m[f'R{GROUP_NAMES[g]}'] = f"{st.loc[grp(g), 'R_all'].mean():.3f}"
    lines = ['% Generated by src/10_paper_tables.py from results/*.json. Do not edit by hand.']
    lines += [f'\\newcommand{{\\{k}}}{{{v}}}' for k, v in m.items()]
    (OUT / 'numbers.tex').write_text('\n'.join(lines) + '\n', encoding='utf-8')

    def row(name, r, p):
        lam = '--' if pd.isna(r.get('lambda', np.nan)) else f"{r['lambda']:.2f}"
        return f"{name} & {num(r['beta'])} & [{num(r['ci95_low'])}, {num(r['ci95_high'])}] & {lam} & {p} \\\\"
    write_macro('PrimaryRows', [
        row('P1: all genes', s['P1'], f"{s['P1']['p_boot']:.3f}"),
        row('P2: GenAge genes', s['P2_selfcontained'], f"{s['P2_selfcontained']['p_boot']:.3f}"),
        f"P2 competitive: GenAge vs 1,000 background sets & -- & -- & -- & {s['P2_competitive']['p_empirical_two_sided']:.3f}$^{{a}}$ \\\\",
        '\\midrule',
        row('Divergence from ancestor per site (positive control)', s['secondary']['R'], f"{s['secondary']['R']['p_boot']:.3f}"),
        row('Net $\\Delta\\Delta G$ per site', s['secondary']['N'], f"{s['secondary']['N']['p_boot']:.3f}")])

    names = {'lambda_fixed_1': '$\\lambda$ fixed at 1', 'variant_no_mask': 'No window mask',
             'variant_plddt90': 'Sites with pLDDT $\\geq$ 90', 'variant_buried': 'Buried sites only',
             'variant_exposed': 'Exposed sites only', 'reptiles_only': 'Reptiles only', 'amphibians_only': 'Amphibians only',
             'without_turtles': 'Turtles excluded', 'squamates_only': 'Squamates only',
             'anage_mass_only': 'Species with AnAge mass only', 'no_mass_covariate': 'No mass covariate',
             'adjust_log10_R': 'Adjusted for log$_{10}$ divergence from ancestor', 'background_genes_only': 'Background genes only',
             'OLS_HC3': 'OLS with HC3 errors (no phylogeny)'}
    rows = []
    for k, r in sens.iterrows():
        if k.startswith('coverage') or pd.isna(r.get('beta')):
            continue
        lam = '--' if pd.isna(r['lambda']) else f"{r['lambda']:.2f}"
        rows.append(f"{names.get(k, k)} & {int(r.n)} & {num(r.beta)} & [{num(r.ci95_low)}, {num(r.ci95_high)}] & {lam} & {pval(r.p_wald)} \\\\")
    write_macro('SensitivityRows', rows)

    sp = pd.read_csv(DATA / 'species_final.csv').set_index('binomial')
    rows = []
    order = {'Testudines': 0, 'Crocodylia': 1, 'Lepidosauria': 2, 'Amphibia': 3}
    st['o'] = st.group.map(order)
    for r in st.sort_values(['o', 'binomial']).itertuples():
        name = sp.loc[r.binomial, 'anage_binomial'].replace('_', ' ').capitalize()
        rows.append(f"\\emph{{{name}}} & {GROUP_NAMES[r.group]} & {r.max_lifespan_yr:g} & {r.mass_g:,.0f} & "
                    f"{r.mass_source} & {r.M_all:.3f} & {r.R_all:.3f} \\\\")
    write_macro('SpeciesRows', rows)
    print('written to', OUT)


if __name__ == '__main__':
    main()
