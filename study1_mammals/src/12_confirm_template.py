"""CONFIRMATORY tests of plan amendment E5 (written and hashed before any platypus template data were retrieved).

Common genes: valid models for both outgroup templates (opossum, platypus). Common sites: pLDDT >= 70 in the human,
opossum and platypus models, aligned in all three, not window-masked. On these identical sites M_s is computed with
the human, opossum and platypus templates.

H5a: Delta_human - Delta_platypus < 0, where Delta_T = (mean M primates - mean M other placentals) / SD(M);
     confirmed if the upper bound of the 95% gene-bootstrap CI is < 0.
H5b: r(M_opossum, M_platypus) - r(M_human, M_platypus) > 0; confirmed if the lower bound of the 95% CI is > 0.
H5c: P1 with the platypus-template metric: replicated null if bootstrap P >= 0.025 and the 95% CI of beta includes 0.
"""
import json

import numpy as np
import pandas as pd
from Bio import Phylo

from common import DATA, RAW, RESULTS
from pgls import PGLS, bootstrap_lrt, summarize, tree_covariance

AA20 = 'ACDEFGHIKLMNPQRSTVWY'
IDX = {a: i for i, a in enumerate(AA20)}
WINDOW, MAX_DIFF = 10, 5
TEMPLATES = {'opossum': 'monodelphis_domestica', 'platypus': 'ornithorhynchus_anatinus'}
SEED_BOOT, SEED_P1, SEED_COMP = 20261059, 20261060, 20261061
B_GENES = 1000


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path, encoding='utf-8'):
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def window_mask(diff):
    mask = np.zeros(len(diff), bool)
    c = np.concatenate([[0], np.cumsum(diff)])
    for w in range(0, max(1, len(diff) - WINDOW + 1)):
        if c[min(w + WINDOW, len(diff))] - c[w] >= MAX_DIFF:
            mask[w:w + WINDOW] = True
    return mask


def gene_sums(g, species):
    z = np.load(DATA / 'energies' / f'{g}.npz')
    P, anc, plddt = z['probs'], str(z['ancestor']), z['plddt']
    E = {'human': z['e_anc']}
    ok = plddt >= 70
    for key, t in TEMPLATES.items():
        o = np.load(DATA / f'energies_{t}' / f'{g}.npz')
        E[key] = o['e']
        ok &= o['template_ok']
    expected = {k: np.nansum(P * np.nan_to_num(v), axis=1) for k, v in E.items()}
    seqs = read_fasta(DATA / 'alignments' / f'{g}.fasta')
    out = {k: np.full(len(species), np.nan) for k in ['human', 'opossum', 'platypus', 'K']}
    for j, sp in enumerate(species):
        s = seqs.get(sp)
        if s is None:
            continue
        code = np.array([IDX.get(c.replace('U', 'C'), -1) for c in s])
        present = code >= 0
        masked = window_mask((present & (np.array(list(s)) != np.array(list(anc)))).astype(int))
        m = present & ok & ~masked
        if not m.any():
            continue
        ci = np.where(present, code, 0)
        rows = np.arange(len(s))
        out['K'][j] = float((1 - P[rows, ci])[m].sum())
        for k in ('human', 'opossum', 'platypus'):
            out[k][j] = float((E[k][rows, ci] - expected[k])[m].sum())
    return out, int(ok.sum())


def offsets(M, prim):
    return (M[prim].mean() - M[~prim].mean()) / M.std(ddof=1)


def main():
    sp = pd.read_csv(DATA / 'species_final.csv')
    sp = sp[sp.eligible & sp.placental].sort_values('ensembl_name').reset_index(drop=True)
    species = sp.ensembl_name.tolist()
    prim = (sp.order == 'Primates').to_numpy()
    sets = [set(f.stem for f in (DATA / f'energies_{t}').glob('*.npz')) for t in TEMPLATES.values()]
    genes = sorted(set.intersection(*sets))
    assert len(genes) >= 500, 'fallback of E5 would apply'
    arrays = {k: [] for k in ['human', 'opossum', 'platypus', 'K']}
    n_sites = 0
    for g in genes:
        o, n = gene_sums(g, species)
        n_sites += n
        for k in arrays:
            arrays[k].append(o[k])
    A = {k: np.column_stack(v) for k, v in arrays.items()}          # species x genes
    M = {k: np.nansum(A[k], 1) / np.nansum(A['K'], 1) for k in ['human', 'opossum', 'platypus']}

    # Gene bootstrap for H5a and H5b
    rng = np.random.default_rng(SEED_BOOT)
    d_off, d_r, boots = [], [], {'off_h': [], 'off_o': [], 'off_p': []}
    for _ in range(B_GENES):
        idx = rng.integers(0, len(genes), len(genes))
        Kb = np.nansum(A['K'][:, idx], 1)
        Mb = {k: np.nansum(A[k][:, idx], 1) / Kb for k in ['human', 'opossum', 'platypus']}
        oh, op, opl = offsets(Mb['human'], prim), offsets(Mb['opossum'], prim), offsets(Mb['platypus'], prim)
        boots['off_h'].append(oh); boots['off_o'].append(op); boots['off_p'].append(opl)
        d_off.append(oh - opl)
        d_r.append(np.corrcoef(Mb['opossum'], Mb['platypus'])[0, 1] - np.corrcoef(Mb['human'], Mb['platypus'])[0, 1])
    ci = lambda v: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))]
    off = {k: float(offsets(M[k], prim)) for k in M}
    r_op_pl = float(np.corrcoef(M['opossum'], M['platypus'])[0, 1])
    r_hu_pl = float(np.corrcoef(M['human'], M['platypus'])[0, 1])
    r_hu_op = float(np.corrcoef(M['human'], M['opossum'])[0, 1])
    h5a_ci, h5b_ci = ci(d_off), ci(d_r)

    # H5c and descriptive P1 slopes on common sites
    tree = Phylo.read(RAW / 'mcc_node.tre', 'nexus')
    model = PGLS(tree_covariance(tree, sp.tree_label.tolist()))
    y = np.log10(sp.max_lifespan_yr.to_numpy())
    lm = np.log10(sp.adult_mass_g.to_numpy())
    zf = lambda x: (x - x.mean()) / x.std(ddof=1)
    design = lambda m: np.column_stack([np.ones(len(y)), lm, zf(m)])
    fit, b = bootstrap_lrt(model, y, design(M['platypus']), 2, 999, SEED_P1)
    h5c = {**summarize(fit, 2), **b}
    slopes = {k: summarize(model.fit(y, design(M[k])), 2) for k in M}

    # Descriptive: GenAge competitive test with the platypus template, on common genes and sites
    gset = pd.read_csv(RESULTS / 'asr_summary.csv').set_index('symbol').gene_set.reindex(genes).to_numpy()
    gen, bg = gset == 'genage', np.where(gset == 'background')[0]
    agg = lambda cols: np.nansum(A['platypus'][:, cols], 1) / np.nansum(A['K'][:, cols], 1)
    beta_gen = model.fit(y, design(agg(np.where(gen)[0])))['beta'][2]
    rng = np.random.default_rng(SEED_COMP)
    null = np.array([model.fit(y, design(agg(rng.choice(bg, gen.sum(), replace=False))))['beta'][2] for _ in range(1000)])

    rank = lambda m: int(pd.Series(m).rank().iloc[species.index('homo_sapiens')])
    out = {
        'amendment': 'E5 (confirmatory; plan hashed before platypus data retrieval)',
        'n_common_genes': len(genes), 'n_common_sites_total': n_sites, 'n_species': len(species),
        'primate_offset_SD': off,
        'H5a': {'statistic': 'Delta_human - Delta_platypus', 'estimate': off['human'] - off['platypus'], 'ci95': h5a_ci,
                'confirmed': h5a_ci[1] < 0},
        'H5b': {'statistic': 'r(opossum,platypus) - r(human,platypus)', 'estimate': r_op_pl - r_hu_pl, 'ci95': h5b_ci,
                'r_opossum_platypus': r_op_pl, 'r_human_platypus': r_hu_pl, 'r_human_opossum': r_hu_op,
                'confirmed': h5b_ci[0] > 0},
        'H5c': {**h5c, 'replicated_null': (h5c['p_boot'] >= 0.025) and (h5c['ci95_low'] <= 0 <= h5c['ci95_high'])},
        'descriptive': {
            'P1_slopes_common_sites_wald': slopes,
            'human_rank_of_57': {k: rank(M[k]) for k in M},
            'genage_competitive_platypus': {'beta_genage': float(beta_gen), 'null_mean': float(null.mean()),
                                            'null_sd': float(null.std()), 'n_genage': int(gen.sum()),
                                            'p_empirical_two_sided': float((np.sum(np.abs(null) >= abs(beta_gen)) + 1) / 1001)},
            'bootstrap_offset_ci95': {k: ci(v) for k, v in boots.items()},
        },
    }
    pd.DataFrame({'ensembl_name': species, 'order': sp.order, 'LQ': sp.LQ, **{f'M_{k}': v for k, v in M.items()}}) \
        .to_csv(RESULTS / 'E5_species_metrics_common_sites.csv', index=False)
    (RESULTS / 'E5_confirmatory.json').write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
