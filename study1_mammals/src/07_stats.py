"""Primary tests P1/P2, secondary and sensitivity analyses (plan section 7)."""
import json

import numpy as np
import pandas as pd
import statsmodels.api as sm
from Bio import Phylo

from common import DATA, RAW, RESULTS, SEED, sha256
from pgls import PGLS, bh, bootstrap_lrt, summarize, tree_covariance

B_PRIMARY = 999
N_NULL_SETS = 1000


def z(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std(ddof=1)


class Data:
    def __init__(self):
        sp = pd.read_csv(DATA / 'species_final.csv')
        self.sp = sp[sp.eligible & sp.placental].sort_values('ensembl_name').reset_index(drop=True)
        met = pd.read_csv(RESULTS / 'species_gene_metrics.csv.gz')
        met = met[met.species.isin(self.sp.ensembl_name)]
        genes = pd.read_csv(RESULTS / 'asr_summary.csv')
        ok = set(pd.read_csv(RESULTS / 'thermompnn_summary.csv').query("status == 'ok'").symbol)
        self.genes = genes[genes.symbol.isin(ok)].set_index('symbol')
        self.symbols = sorted(self.genes.index)
        self.mat = {}
        for v in ['main', 'no_mask', 'plddt90', 'buried', 'exposed']:
            for q in 'DKL':
                self.mat[q, v] = (met.pivot(index='species', columns='symbol', values=f'{q}_{v}')
                                  .reindex(index=self.sp.ensembl_name, columns=self.symbols).to_numpy())
        self.gene_set = self.genes.gene_set.reindex(self.symbols).to_numpy()
        self.old = self.genes.old_family.reindex(self.symbols).fillna(False).astype(bool).to_numpy()

    def aggregate(self, cols, variant='main', kind='M'):
        D = np.nansum(self.mat['D', variant][:, cols], axis=1)
        K = np.nansum(self.mat['K', variant][:, cols], axis=1)
        L = np.nansum(self.mat['L', variant][:, cols], axis=1)
        return {'M': D / K, 'N': D / L, 'R': K / L}[kind]


def design(sp, metric, extra=()):
    cols = [np.ones(len(sp)), np.log10(sp.adult_mass_g.to_numpy())] + [z(e) for e in extra] + [z(metric)]
    return np.log10(sp.max_lifespan_yr.to_numpy()), np.column_stack(cols)


def fit_one(tree, sp, metric, extra=(), lam=None, boot=0, seed=SEED):
    y, X = design(sp, metric, extra)
    model = PGLS(tree_covariance(tree, sp.tree_label.tolist()))
    col = X.shape[1] - 1
    if boot:
        fit, b = bootstrap_lrt(model, y, X, col, boot, seed, lam)
        return {**summarize(fit, col), **b, 'n': len(y)}
    return {**summarize(model.fit(y, X, lam), col), 'n': len(y)}


def main():
    d = Data()
    sp = d.sp
    tree = Phylo.read(RAW / 'mcc_node.tre', 'nexus')
    allg = np.ones(len(d.symbols), bool)
    gen = d.gene_set == 'genage'
    bg = d.gene_set == 'background'
    out = {'n_species': len(sp), 'n_genes': int(allg.sum()), 'n_genage': int(gen.sum()), 'n_background': int(bg.sum()),
           'inputs_sha256': {'species_gene_metrics': sha256(RESULTS / 'species_gene_metrics.csv.gz'),
                             'plan': sha256(DATA.parent / 'PLAN_ANALISIS.md')}}
    species_table = sp[['ensembl_name', 'tree_label', 'order', 'max_lifespan_yr', 'adult_mass_g', 'LQ', 'data_quality']].copy()
    for name, cols in [('all', allg), ('genage', gen), ('background', bg)]:
        for kind in 'MNR':
            species_table[f'{kind}_{name}'] = d.aggregate(cols, kind=kind)
    species_table['LQ_group'] = pd.qcut(species_table.LQ, 3, labels=['corta', 'intermedia', 'larga'])
    species_table.to_csv(RESULTS / 'species_metrics.csv', index=False)

    # P1
    out['P1'] = fit_one(tree, sp, species_table.M_all, boot=B_PRIMARY, seed=SEED + 1)
    # P2 self-contained and competitive
    out['P2_selfcontained'] = fit_one(tree, sp, species_table.M_genage, boot=B_PRIMARY, seed=SEED + 2)
    beta_obs = out['P2_selfcontained']['beta']
    rng = np.random.default_rng(SEED + 3)
    bg_idx = np.where(bg)[0]
    y, _ = design(sp, species_table.M_all)
    model = PGLS(tree_covariance(tree, sp.tree_label.tolist()))
    null_betas = []
    for _ in range(N_NULL_SETS):
        cols = np.zeros(len(d.symbols), bool)
        cols[rng.choice(bg_idx, gen.sum(), replace=False)] = True
        _, X = design(sp, d.aggregate(cols))
        null_betas.append(model.fit(y, X)['beta'][2])
    null_betas = np.array(null_betas)
    np.save(RESULTS / 'P2_null_betas.npy', null_betas)
    out['P2_competitive'] = {'beta_genage': beta_obs, 'null_mean': float(null_betas.mean()), 'null_sd': float(null_betas.std()),
                             'p_empirical_two_sided': float((np.sum(np.abs(null_betas) >= abs(beta_obs)) + 1) / (N_NULL_SETS + 1)),
                             'set_size': int(gen.sum()), 'n_null_sets': N_NULL_SETS}
    out['primary_alpha_bonferroni'] = 0.025

    # Secondary
    sec = {}
    for kind in 'NR':
        sec[f'{kind}_all'] = fit_one(tree, sp, species_table[f'{kind}_all'], boot=B_PRIMARY, seed=SEED + 10 + ord(kind))
    mm = PGLS(tree_covariance(tree, sp.tree_label.tolist()))
    Xm = np.column_stack([np.ones(len(sp)), np.log10(sp.adult_mass_g)])
    sec['nearly_neutral_M_on_log10mass'] = summarize(mm.fit(z(species_table.M_all), Xm), 1)
    sec['nearly_neutral_M_on_log10lifespan_and_mass'] = summarize(
        mm.fit(z(species_table.M_all), np.column_stack([Xm, np.log10(sp.max_lifespan_yr)])), 2)
    out['secondary'] = sec

    # Sensitivity of P1 (analytic p values; diagnostics, not new discoveries)
    sens = []
    def add(name, sub, metric, extra=(), lam=None, tr=tree):
        try:
            r = fit_one(tr, sub, metric, extra, lam)
        except ValueError as e:
            r = {'error': str(e)}
        sens.append({'analysis': name, **r})
    add('lambda_fixed_1', sp, species_table.M_all, lam=1.0)
    for v in ['no_mask', 'plddt90', 'buried', 'exposed']:
        add(f'variant_{v}', sp, d.aggregate(allg, v))
    for label, keep in [('exclude_primates', sp.order != 'Primates'), ('exclude_human', sp.ensembl_name != 'homo_sapiens'),
                        ('exclude_bats', sp.order != 'Chiroptera'), ('high_quality_only', sp.data_quality == 'high'),
                        ('exclude_rodents', sp.order != 'Rodentia')]:
        sub = sp[keep.to_numpy()].reset_index(drop=True)
        add(label, sub, species_table.M_all[keep.to_numpy()].to_numpy())
    add('adjust_log10_R', sp, species_table.M_all, extra=[np.log10(species_table.R_all)])
    add('background_genes_only', sp, species_table.M_background)
    yq = np.log10(sp.LQ.to_numpy())
    fq = mm.fit(yq, np.column_stack([np.ones(len(sp)), z(species_table.M_all)]))
    sens.append({'analysis': 'response_log10_LQ_no_mass', **summarize(fq, 1), 'n': len(sp)})
    yv, Xv = design(sp, species_table.M_all)
    ols = sm.OLS(yv, Xv).fit(cov_type='HC3')
    sens.append({'analysis': 'OLS_HC3', 'beta': float(ols.params[-1]), 'ci95_low': float(ols.conf_int()[-1][0]),
                 'ci95_high': float(ols.conf_int()[-1][1]), 'p_wald': float(ols.pvalues[-1]), 'n': len(yv)})
    keep_idx = set(np.round(np.linspace(0, 99, 20)).astype(int))
    for i, t in enumerate(Phylo.parse(RAW / 'posterior100.trees', 'nexus')):
        if i in keep_idx:
            add(f'posterior_tree_{i}', sp, species_table.M_all, tr=t)
        if i >= 99:
            break
    sens = pd.DataFrame(sens)
    sens.to_csv(RESULTS / 'sensitivity_P1.csv', index=False)

    # Exploratory per-gene tests
    rows = []
    D, K = d.mat['D', 'main'], d.mat['K', 'main']
    for j, g in enumerate(d.symbols):
        m = D[:, j] / K[:, j]
        ok = np.isfinite(m) & (K[:, j] >= 0.5)
        if ok.sum() < 30 or np.nanstd(m[ok]) == 0:
            continue
        sub = sp[ok].reset_index(drop=True)
        try:
            r = fit_one(tree, sub, m[ok])
        except (ValueError, np.linalg.LinAlgError):
            continue
        rows.append({'symbol': g, 'gene_set': d.gene_set[j], 'old_family': bool(d.old[j]), **r})
    genes = pd.DataFrame(rows)
    genes['q_bh'] = bh(genes.p_wald)
    genes.to_csv(RESULTS / 'per_gene_exploratory.csv', index=False)
    out['per_gene'] = {'tested': len(genes), 'q_lt_0.05': int((genes.q_bh < 0.05).sum()),
                       'min_q': float(genes.q_bh.min()), 'fraction_beta_positive': float((genes.beta > 0).mean())}
    se = out['P1']['se']
    out['P1_detectable_beta_80pct_power'] = float((1.96 + 0.8416) * se)
    out['P1_detectable_lifespan_percent_per_SD'] = float(100 * (10 ** ((1.96 + 0.8416) * se) - 1))
    (RESULTS / 'stats_summary.json').write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps({k: out[k] for k in ['n_species', 'n_genes', 'n_genage', 'P1', 'P2_selfcontained', 'P2_competitive']},
                     indent=1, default=float))


if __name__ == '__main__':
    main()
