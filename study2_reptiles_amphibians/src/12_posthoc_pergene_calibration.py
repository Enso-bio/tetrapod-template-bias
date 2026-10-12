"""POST HOC amendment E2(b): calibration of the exploratory per-gene tests by simulation under the null.

The null model (log10 lifespan ~ log10 mass, lambda by ML) is fitted to all species and 1,000 responses are simulated
from it. For each simulated response the per-gene tests of 07_stats.py are repeated exactly (same species per gene,
lambda by ML, Wald test); a variant fixes lambda at the null-model value.
"""
import importlib
import json
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd
from Bio import Phylo

from common import DATA, RESULTS
from pgls import PGLS, bh, tree_covariance

stats = importlib.import_module('07_stats')
B, SEED = 1000, 20261082
_G = {}


def build_tests():
    """Per-gene designs and PGLS models, with the same inclusion rules as the per-gene tests of 07_stats.py."""
    d = stats.Data()
    sp = d.sp
    tree = Phylo.read(DATA / 'timetree_upgma.nwk', 'newick')
    lm = np.log10(sp.mass_g.to_numpy())
    D, K = d.mat['D', 'main'], d.mat['K', 'main']
    models, tests = {}, []
    for j, g in enumerate(d.ids):
        m = D[:, j] / K[:, j]
        ok = np.isfinite(m) & (K[:, j] >= 0.5)
        if ok.sum() < 25 or np.nanstd(m[ok]) == 0:
            continue
        key = ok.tobytes()
        if key not in models:
            models[key] = PGLS(tree_covariance(tree, sp.binomial[ok].tolist()))
        X = np.column_stack([np.ones(ok.sum()), lm[ok], stats.z(m[ok])])
        tests.append((g, ok, key, models[key].project(X)))
    return sp, tree, models, tests


def _init(models, tests, lam0):
    _G.update(models=models, tests=tests, lam0=lam0)


def _run(ys):
    """Wald P of every per-gene test for each simulated response (lambda by ML and fixed at lam0)."""
    models, tests, lam0 = _G['models'], _G['tests'], _G['lam0']
    p_ml = np.empty((len(ys), len(tests)), np.float32)
    p_fx = np.empty_like(p_ml)
    for b, y in enumerate(ys):
        for k, (_, ok, key, Xp) in enumerate(tests):
            yp = models[key].project(y[ok])
            p_ml[b, k] = models[key].fit_projected(yp, Xp)['p'][2]
            p_fx[b, k] = models[key].fit_projected(yp, Xp, lam0)['p'][2]
    return p_ml, p_fx


def main():
    sp, tree, models, tests = build_tests()
    y = np.log10(sp.max_lifespan_yr.to_numpy())
    X0 = np.column_stack([np.ones(len(sp)), np.log10(sp.mass_g.to_numpy())])
    C = tree_covariance(tree, sp.binomial.tolist())
    null = PGLS(C).fit(y, X0)
    lam0, s2 = null['lambda'], null['sigma2_ml']
    V = lam0 * C + (1 - lam0) * np.diag(np.diag(C))
    L = np.linalg.cholesky(V)
    rng = np.random.default_rng(SEED)
    Y = X0 @ null['beta'] + np.sqrt(s2) * rng.standard_normal((B, len(sp))) @ L.T

    _init(models, tests, lam0)
    obs_ml, obs_fx = (p[0] for p in _run([y]))
    prev = pd.read_csv(RESULTS / 'per_gene_exploratory.csv').set_index('human_gene').p_wald
    genes = [g for g, *_ in tests]
    assert np.allclose(prev.loc[genes].to_numpy(), obs_ml, rtol=1e-4, atol=1e-6), 'observed per-gene P not reproduced'

    chunks = np.array_split(Y, 40)
    with ProcessPoolExecutor(max_workers=18, initializer=_init, initargs=(models, tests, lam0)) as pool:
        parts = list(pool.map(_run, chunks))
    P_ml = np.vstack([p for p, _ in parts])
    P_fx = np.vstack([p for _, p in parts])
    np.savez_compressed(RESULTS / 'posthoc_E2b_null_pvalues.npz', p_ml=P_ml, p_fixed=P_fx, genes=np.array(genes, dtype=object))

    def summary(obs, sim):
        n_obs = int((obs < 0.05).sum())
        counts = (sim < 0.05).sum(axis=1)
        emp = (1 + (sim <= obs[None, :]).sum(axis=0)) / (B + 1)
        return {'observed_nominal': n_obs,
                'null_mean': float(counts.mean()), 'null_sd': float(counts.std(ddof=1)),
                'null_interval95': [float(np.percentile(counts, 2.5)), float(np.percentile(counts, 97.5))],
                'null_false_positive_rate': float(counts.mean() / len(obs)),
                'p_calibrated_count': float((1 + (counts >= n_obs).sum()) / (B + 1)),
                'genes_empirical_p_lt_0.05': int((emp < 0.05).sum()),
                'genes_empirical_q_lt_0.05': int((bh(emp) < 0.05).sum()),
                'min_empirical_p': float(emp.min())}

    per_gene = pd.read_csv(RESULTS / 'per_gene_exploratory.csv')
    out = {'amendment': 'E2(b) (post hoc, exploratory)', 'n_tests': len(tests), 'B': B,
           'null_model': {'lambda': float(lam0), 'sigma2_ml': float(s2), 'beta': null['beta'].tolist()},
           'observed_fraction_lambda_lt_0.05': float((per_gene['lambda'] < 0.05).mean()),
           'lambda_ML': summary(obs_ml, P_ml), 'lambda_fixed_null': summary(obs_fx, P_fx)}
    (RESULTS / 'posthoc_E2b_pergene_calibration.json').write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
