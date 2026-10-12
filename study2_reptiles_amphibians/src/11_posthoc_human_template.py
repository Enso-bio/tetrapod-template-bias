"""POST HOC amendment E2(a): follow-up of the nominal association with the human template in reptiles.

Uses the species metrics of the template sensitivity analysis (300 genes, sites shared by the zebrafish, coelacanth
and human templates; 08_template_sensitivity.py) for the 29 reptiles.
"""
import json

import numpy as np
import pandas as pd
from Bio import Phylo

from common import DATA, RESULTS
from pgls import PGLS, bootstrap_lrt, summarize, tree_covariance

REPTILES = {'Testudines', 'Crocodylia', 'Lepidosauria'}
ARCHELOSAURS = {'Testudines', 'Crocodylia'}
SEEDS = {'homo_sapiens': 20261080, 'zebrafish': 20261081}


def z(x):
    x = np.asarray(x, float)
    return (x - x.mean()) / x.std(ddof=1)


def fit(tree, sub, metric, extra=(), boot=0, seed=0):
    y = np.log10(sub.max_lifespan_yr.to_numpy())
    X = np.column_stack([np.ones(len(sub)), np.log10(sub.mass_g.to_numpy())] + list(extra) + [z(metric)])
    model = PGLS(tree_covariance(tree, sub.binomial.tolist()))
    col = X.shape[1] - 1
    if boot:
        f, b = bootstrap_lrt(model, y, X, col, boot, seed)
        return {**summarize(f, col), **b, 'n': len(y)}
    return {**summarize(model.fit(y, X), col), 'n': len(y)}


def main():
    sp = pd.read_csv(DATA / 'species_final.csv')
    sp = sp[sp.eligible].sort_values('binomial').reset_index(drop=True)
    M = pd.read_csv(RESULTS / 'template_sensitivity_species_M.csv', index_col=0).loc[sp.binomial]
    R = pd.read_csv(RESULTS / 'species_metrics.csv').set_index('binomial').R_all.loc[sp.binomial]
    rept = sp.group.isin(REPTILES).to_numpy()
    sub = sp[rept].reset_index(drop=True)
    Mr, Rr = M[rept], np.log10(R[rept].to_numpy())
    arch = sub.group.isin(ARCHELOSAURS).to_numpy().astype(float)
    tree = Phylo.read(DATA / 'timetree_upgma.nwk', 'newick')
    out = {'amendment': 'E2(a) (post hoc, exploratory)', 'n_reptiles': len(sub),
           'n_by_clade': sub.group.value_counts().to_dict(), 'templates': {}}
    for t in ['homo_sapiens', 'zebrafish']:
        m = Mr[t].to_numpy()
        res = {'original_wald': fit(tree, sub, m),
               'bootstrap_lrt': fit(tree, sub, m, boot=999, seed=SEEDS[t]),
               'clade_covariate': fit(tree, sub, m, extra=[arch]),
               'adjusted_log10_R': fit(tree, sub, m, extra=[z(Rr)])}
        for label, keep in [('squamates_only', arch == 0), ('turtles_crocodilians_only', arch == 1)]:
            s2 = sub[keep].reset_index(drop=True)
            res[label] = fit(tree, s2, m[keep])
        zm = z(m)
        res['clade_mean_z_metric'] = {'turtles_crocodilians': float(zm[arch == 1].mean()), 'squamates': float(zm[arch == 0].mean())}
        res['r_metric_log10_R'] = float(np.corrcoef(m, Rr)[0, 1])
        out['templates'][t] = res
    diff = z(Mr.homo_sapiens) - z(Mr.zebrafish)
    out['human_minus_zebrafish_z'] = {
        'r_with_log10_R': float(np.corrcoef(diff, Rr)[0, 1]),
        'mean_turtles_crocodilians': float(diff[arch == 1].mean()), 'mean_squamates': float(diff[arch == 0].mean())}
    out['log10_R_by_clade'] = {'turtles_crocodilians': float(Rr[arch == 1].mean()), 'squamates': float(Rr[arch == 0].mean())}
    (RESULTS / 'posthoc_E2a_human_template.json').write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
