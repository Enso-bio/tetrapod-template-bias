"""Phylogenetic GLS with Pagel's lambda (ML), parametric-bootstrap likelihood-ratio test and tree covariance."""
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import t as tdist


def tree_covariance(tree, labels):
    """C[i,j] = shared path length from the MRCA of `labels` to MRCA(i,j); rows ordered as labels.

    Starting at the MRCA of the analysed species drops any stem/outgroup branch shared by all of them.
    """
    tips = {t.name: t for t in tree.get_terminals()}
    missing = set(labels) - set(tips)
    if missing:
        raise ValueError(f'labels absent from tree: {sorted(missing)}')
    start = tree.common_ancestor([tips[l] for l in labels])
    paths = {lab: start.get_path(tips[lab]) for lab in labels}
    n = len(labels)
    C = np.zeros((n, n))
    for i, a in enumerate(labels):
        for j in range(i, n):
            b = labels[j]
            shared = 0.0
            for x, y in zip(paths[a], paths[b]):
                if x is not y:
                    break
                shared += float(x.branch_length or 0)
            C[i, j] = C[j, i] = shared
    if np.min(np.diag(C)) <= 0:
        raise ValueError('zero root-to-tip distance')
    return C / np.mean(np.diag(C))


class PGLS:
    """V(lambda) = lambda*C + (1-lambda)*diag(C) = D U [lambda*E + (1-lambda) I] U' D, with R = D^-1 C D^-1 = U E U'."""

    def __init__(self, C):
        self.n = len(C)
        self.d = np.sqrt(np.diag(C))
        e, U = np.linalg.eigh(C / np.outer(self.d, self.d))
        if e.min() < -1e-8:
            raise ValueError('covariance not positive semidefinite')
        self.e, self.U = np.maximum(e, 1e-12), U
        self.logdet_d = 2 * np.sum(np.log(self.d))

    def project(self, a):
        a = np.asarray(a, float)
        return self.U.T @ (a / self.d if a.ndim == 1 else a / self.d[:, None])

    def _calc(self, yp, Xp, lam, detail=False):
        n, p = Xp.shape
        var = (1 - lam) + lam * self.e
        w = 1 / np.sqrt(var)
        xw, yw = Xp * w[:, None], yp * w
        beta, *_ = np.linalg.lstsq(xw, yw, rcond=None)
        r = yw - xw @ beta
        sse = float(r @ r)
        obj = n * np.log(sse / n) + np.sum(np.log(var)) + self.logdet_d      # -2 logL up to a constant
        if not detail:
            return obj
        cov = np.linalg.inv(xw.T @ xw) * sse / (n - p)
        se = np.sqrt(np.diag(cov))
        tval = beta / se
        return {'lambda': float(lam), 'objective': float(obj), 'beta': beta, 'se': se, 'df': n - p,
                'p': 2 * tdist.sf(np.abs(tval), n - p), 'sigma2_ml': sse / n, 'var_eigen': var}

    def fit_projected(self, yp, Xp, lam=None):
        if lam is None:
            grid = np.linspace(0, 1, 41)
            objs = [self._calc(yp, Xp, g) for g in grid]
            k = int(np.argmin(objs))
            lo, hi = grid[max(k - 1, 0)], grid[min(k + 1, len(grid) - 1)]
            opt = minimize_scalar(lambda v: self._calc(yp, Xp, v), bounds=(lo, hi), method='bounded',
                                  options={'xatol': 1e-4})
            cands = [grid[k], float(opt.x)]
            lam = min(cands, key=lambda v: self._calc(yp, Xp, v))
        return self._calc(yp, Xp, lam, detail=True)

    def fit(self, y, X, lam=None):
        return self.fit_projected(self.project(y), self.project(X), lam)


def summarize(fit, col):
    b, se = fit['beta'][col], fit['se'][col]
    c = tdist.ppf(0.975, fit['df']) * se
    return {'beta': float(b), 'se': float(se), 'ci95_low': float(b - c), 'ci95_high': float(b + c),
            'p_wald': float(fit['p'][col]), 'lambda': fit['lambda']}


def bootstrap_lrt(model, y, X, col, B=999, seed=0, lam_fixed=None):
    """Parametric bootstrap under the null (column `col` removed); lambda re-estimated in each replicate."""
    yp, Xp = model.project(y), model.project(X)
    X0 = np.delete(Xp, col, axis=1)
    null, full = model.fit_projected(yp, X0, lam_fixed), model.fit_projected(yp, Xp, lam_fixed)
    observed = max(0.0, null['objective'] - full['objective'])
    rng = np.random.default_rng(seed)
    mean = X0 @ null['beta']
    sd = np.sqrt(null['sigma2_ml'] * null['var_eigen'])
    exceed = 0
    for _ in range(B):
        sim = mean + sd * rng.standard_normal(len(y))
        lr = max(0.0, model.fit_projected(sim, X0, lam_fixed)['objective'] - model.fit_projected(sim, Xp, lam_fixed)['objective'])
        exceed += lr >= observed - 1e-10
    p = (exceed + 1) / (B + 1)
    return full, {'LRT': observed, 'p_boot': p, 'boot_exceed': exceed, 'boot_B': B, 'null_lambda': null['lambda'],
                  'p_boot_mc_se': float(np.sqrt(p * (1 - p) / (B + 1)))}


def bh(p):
    p = np.asarray(p, float)
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[order] = np.minimum(q, 1)
    return out


def self_test():
    """Eigen-transformed PGLS equals direct GLS; lambda=0 equals weighted OLS; toy tree covariance."""
    from Bio import Phylo
    import io
    rng = np.random.default_rng(1)
    A = rng.normal(size=(25, 25))
    C = A @ A.T + 25 * np.eye(25)
    y, X = rng.normal(size=25), np.column_stack([np.ones(25), rng.normal(size=(25, 2))])
    m = PGLS(C)
    for lam in (0.0, 0.37, 1.0):
        f = m.fit(y, X, lam)
        V = lam * C + (1 - lam) * np.diag(np.diag(C))
        Vi = np.linalg.inv(V)
        b = np.linalg.solve(X.T @ Vi @ X, X.T @ Vi @ y)
        r = y - X @ b
        obj = 25 * np.log(r @ Vi @ r / 25) + np.linalg.slogdet(V)[1]
        assert np.allclose(b, f['beta'], atol=1e-9) and np.isclose(obj, f['objective'], atol=1e-8)
    t = Phylo.read(io.StringIO('((a:1,b:1):2,c:3);'), 'newick')
    assert np.allclose(tree_covariance(t, ['a', 'b', 'c']), np.array([[3, 2, 0], [2, 3, 0], [0, 0, 3]]) / 3)
    assert np.allclose(bh([0.01, 0.04, 0.03, 1]), [0.04, 0.04 * 4 / 3, 0.04 * 4 / 3, 1])
    grid = np.linspace(0, 1, 2001)
    f = m.fit(y, X)
    assert f['objective'] <= min(m.fit(y, X, g)['objective'] for g in grid) + 1e-6
    return True


if __name__ == '__main__':
    print('self-test passed:', self_test())
