"""Primary tests P1/P2, secondary and sensitivity analyses (plan section 6)."""
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from Bio import Phylo

from common import DATA, RESULTS, SEED, sha256
from pgls import PGLS, bh, bootstrap_lrt, summarize, tree_covariance

B_PRIMARY = 999
N_NULL_SETS = 1000
REPTILES = {'Testudines', 'Crocodylia', 'Lepidosauria'}


def z(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std(ddof=1)


class Data:
    def __init__(self):
        sp = pd.read_csv(DATA / 'species_final.csv')
        self.sp = sp[sp.eligible].sort_values('binomial').reset_index(drop=True)
        met = pd.read_csv(RESULTS / 'species_gene_metrics.csv.gz')
        met = met[met.species.isin(self.sp.binomial)]
        genes = pd.read_csv(RESULTS / 'asr_summary.csv')
        ok = set(pd.read_csv(RESULTS / 'thermompnn_summary.csv').query("status == 'ok'").human_gene)
        self.genes = genes[genes.human_gene.isin(ok)].set_index('human_gene')
        self.ids = sorted(self.genes.index)
        self.mat = {}
        for v in ['main', 'no_mask', 'plddt90', 'buried', 'exposed']:
            for q in 'DKL':
                self.mat[q, v] = (met.pivot(index='species', columns='human_gene', values=f'{q}_{v}')
                                  .reindex(index=self.sp.binomial, columns=self.ids).to_numpy())
        self.gene_set = self.genes.gene_set.reindex(self.ids).to_numpy()

    def aggregate(self, cols, variant='main', kind='M'):
        D = np.nansum(self.mat['D', variant][:, cols], axis=1)
        K = np.nansum(self.mat['K', variant][:, cols], axis=1)
        L = np.nansum(self.mat['L', variant][:, cols], axis=1)
        return {'M': D / K, 'N': D / L, 'R': K / L}[kind]


def design(sp, metric, extra=(), mass=True):
    cols = [np.ones(len(sp))] + ([np.log10(sp.mass_g.to_numpy())] if mass else []) + [z(e) for e in extra] + [z(metric)]
    return np.log10(sp.max_lifespan_yr.to_numpy()), np.column_stack(cols)


def fit_one(tree, sp, metric, extra=(), lam=None, boot=0, seed=SEED, mass=True):
    y, X = design(sp, metric, extra, mass)
    model = PGLS(tree_covariance(tree, sp.binomial.tolist()))
    col = X.shape[1] - 1
    if boot:
        fit, b = bootstrap_lrt(model, y, X, col, boot, seed, lam)
        return {**summarize(fit, col), **b, 'n': len(y)}
    return {**summarize(model.fit(y, X, lam), col), 'n': len(y)}


def main():
    d = Data()
    sp = d.sp
    tree = Phylo.read(DATA / 'timetree_upgma.nwk', 'newick')
    allg = np.ones(len(d.ids), bool)
    gen, bg = d.gene_set == 'genage', d.gene_set == 'background'
    out = {'n_species': len(sp), 'n_by_group': sp.group.value_counts().to_dict(), 'n_genes': int(allg.sum()),
           'n_genage': int(gen.sum()), 'n_background': int(bg.sum()),
           'inputs_sha256': {'species_gene_metrics': sha256(RESULTS / 'species_gene_metrics.csv.gz'),
                             'plan': sha256(DATA.parent / 'PLAN_ANALISIS.md')}}
    st = sp[['binomial', 'organism', 'group', 'order', 'max_lifespan_yr', 'mass_g', 'mass_source', 'data_quality',
             'body_temp_c', 'range_mean_temp_c']].copy()
    for name, cols in [('all', allg), ('genage', gen), ('background', bg)]:
        for kind in 'MNR':
            st[f'{kind}_{name}'] = d.aggregate(cols, kind=kind)
    st.to_csv(RESULTS / 'species_metrics.csv', index=False)

    out['P1'] = fit_one(tree, sp, st.M_all, boot=B_PRIMARY, seed=SEED + 1)
    out['P2_selfcontained'] = fit_one(tree, sp, st.M_genage, boot=B_PRIMARY, seed=SEED + 2)
    rng = np.random.default_rng(SEED + 3)
    y, _ = design(sp, st.M_all)
    model = PGLS(tree_covariance(tree, sp.binomial.tolist()))
    null = []
    for _ in range(N_NULL_SETS):
        cols = np.zeros(len(d.ids), bool)
        cols[rng.choice(np.where(bg)[0], gen.sum(), replace=False)] = True
        null.append(model.fit(y, design(sp, d.aggregate(cols))[1])['beta'][2])
    null = np.array(null)
    np.save(RESULTS / 'P2_null_betas.npy', null)
    b0 = out['P2_selfcontained']['beta']
    out['P2_competitive'] = {'beta_genage': b0, 'null_mean': float(null.mean()), 'null_sd': float(null.std()),
                             'p_empirical_two_sided': float((np.sum(np.abs(null) >= abs(b0)) + 1) / (N_NULL_SETS + 1)),
                             'set_size': int(gen.sum()), 'n_null_sets': N_NULL_SETS}
    out['primary_alpha_bonferroni'] = 0.025

    sec = {k: fit_one(tree, sp, st[f'{k}_all'], boot=B_PRIMARY, seed=SEED + 10 + i) for i, k in enumerate('NR')}
    mm = PGLS(tree_covariance(tree, sp.binomial.tolist()))
    sec['nearly_neutral_M_on_log10mass'] = summarize(
        mm.fit(z(st.M_all), np.column_stack([np.ones(len(sp)), np.log10(sp.mass_g)])), 1)
    out['secondary'] = sec

    sens = []
    def add(name, keep=None, metric=None, extra=(), lam=None, mass=True):
        keep = np.ones(len(sp), bool) if keep is None else np.asarray(keep)
        m = (st.M_all if metric is None else metric)[keep]
        sub = sp[keep].reset_index(drop=True)
        ex = [np.asarray(e)[keep] for e in extra]
        try:
            r = fit_one(tree, sub, np.asarray(m), ex, lam, mass=mass)
        except (ValueError, np.linalg.LinAlgError) as e:
            r = {'error': str(e)}
        sens.append({'analysis': name, **r})
    add('lambda_fixed_1', lam=1.0)
    for v in ['no_mask', 'plddt90', 'buried', 'exposed']:
        add(f'variant_{v}', metric=pd.Series(d.aggregate(allg, v)))
    rept = sp.group.isin(REPTILES).to_numpy()
    add('reptiles_only', rept)
    add('amphibians_only', (sp.group == 'Amphibia').to_numpy())
    add('without_turtles', (sp.group != 'Testudines').to_numpy())
    add('squamates_only', (sp.order == 'Squamata').to_numpy())
    add('anage_mass_only', (sp.mass_source == 'AnAge').to_numpy())
    add('no_mass_covariate', mass=False)
    add('adjust_log10_R', extra=[np.log10(st.R_all)])
    add('background_genes_only', metric=st.M_background)
    for col in ['body_temp_c', 'range_mean_temp_c']:
        cover = st.loc[rept, col].notna().mean()
        keep = rept & st[col].notna().to_numpy()
        if cover >= 0.70:
            add(f'reptiles_adjust_{col}', keep, extra=[st[col].to_numpy()])
        sens.append({'analysis': f'coverage_{col}_reptiles', 'coverage': float(cover)})
    yv, Xv = design(sp, st.M_all)
    ols = sm.OLS(yv, Xv).fit(cov_type='HC3')
    sens.append({'analysis': 'OLS_HC3', 'beta': float(ols.params[-1]), 'ci95_low': float(ols.conf_int()[-1][0]),
                 'ci95_high': float(ols.conf_int()[-1][1]), 'p_wald': float(ols.pvalues[-1]), 'n': len(yv)})
    pd.DataFrame(sens).to_csv(RESULTS / 'sensitivity_P1.csv', index=False)

    rows = []
    D, K = d.mat['D', 'main'], d.mat['K', 'main']
    for j, g in enumerate(d.ids):
        m = D[:, j] / K[:, j]
        ok = np.isfinite(m) & (K[:, j] >= 0.5)
        if ok.sum() < 25 or np.nanstd(m[ok]) == 0:
            continue
        try:
            r = fit_one(tree, sp[ok].reset_index(drop=True), m[ok])
        except (ValueError, np.linalg.LinAlgError):
            continue
        rows.append({'human_gene': g, 'symbol': d.genes.loc[g, 'symbol'], 'gene_set': d.gene_set[j], **r})
    genes = pd.DataFrame(rows)
    genes['q_bh'] = bh(genes.p_wald)
    genes.to_csv(RESULTS / 'per_gene_exploratory.csv', index=False)
    out['per_gene'] = {'tested': len(genes), 'q_lt_0.05': int((genes.q_bh < 0.05).sum()), 'min_q': float(genes.q_bh.min()),
                       'nominal_p_lt_0.05': int((genes.p_wald < 0.05).sum()),
                       'fraction_beta_positive': float((genes.beta > 0).mean())}
    se = out['P1']['se']
    out['P1_detectable_lifespan_percent_per_SD'] = float(100 * (10 ** ((1.96 + 0.8416) * se) - 1))
    (RESULTS / 'stats_summary.json').write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps({k: out[k] for k in ['n_species', 'n_by_group', 'n_genes', 'n_genage', 'P1', 'P2_selfcontained',
                                          'P2_competitive']}, indent=1, default=float))


if __name__ == '__main__':
    main()
