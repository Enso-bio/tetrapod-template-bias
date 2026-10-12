"""Species x gene stability metrics (plan section 6) and predictor sanity controls (C1-C3).

For species s, gene g and eligible site i with residue b:
  ddG_i(s) = E[i,b] - sum_a P(anc_i=a) E[i,a]      (ancestral context; positive = destabilizing)
  K_i(s)   = 1 - P(anc_i=b)
D = sum ddG, K = sum K_i, L = number of sites. Variants: main (pLDDT>=70, window mask), no_mask, plddt90,
buried (RSA<0.25), exposed (RSA>=0.25).
"""
from concurrent.futures import ProcessPoolExecutor

import numpy as np
import pandas as pd

from common import DATA, RESULTS

AA20 = 'ACDEFGHIKLMNPQRSTVWY'
IDX = {a: i for i, a in enumerate(AA20)}
WINDOW, MAX_DIFF = 10, 5


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


def gene_metrics(args):
    symbol, species = args
    z = np.load(DATA / 'energies' / f'{symbol}.npz')
    E, P = z['e_anc'], z['probs']
    anc, plddt, rsa, ss = str(z['ancestor']), z['plddt'], z['rsa'], z['ss']
    expected_anc = (P * E).sum(axis=1)                       # sum_a P(a) E[i,a]
    seqs = read_fasta(DATA / 'alignments' / f'{symbol}.fasta')
    rows = []
    for sp in species:
        s = seqs.get(sp)
        if s is None:
            continue
        code = np.array([IDX.get(c.replace('U', 'C'), -1) for c in s])
        present = code >= 0
        diff = present & (np.array(list(s)) != np.array(list(anc)))
        masked = window_mask(diff.astype(int))
        ci = np.where(present, code, 0)
        ddg = E[np.arange(len(s)), ci] - expected_anc
        k = 1 - P[np.arange(len(s)), ci]
        base = present & (plddt >= 70)
        variants = {'main': base & ~masked, 'no_mask': base, 'plddt90': present & (plddt >= 90) & ~masked,
                    'buried': base & ~masked & (rsa < 0.25), 'exposed': base & ~masked & (rsa >= 0.25)}
        out = {'symbol': symbol, 'species': sp, 'eligible_sites': int((plddt >= 70).sum()),
               'masked_sites': int((masked & base).sum())}
        for name, m in variants.items():
            out[f'D_{name}'] = float(ddg[m].sum())
            out[f'K_{name}'] = float(k[m].sum())
            out[f'L_{name}'] = int(m.sum())
        rows.append(out)

    # Controls on this gene (ancestral MAP context): observed substitutions vs all 19 alternatives.
    a_idx = np.array([IDX[c] for c in anc])
    dd_all = E - E[np.arange(len(anc)), a_idx][:, None]       # ddG(anc_MAP -> x), shape (L, 20)
    alt_mean = (dd_all.sum(axis=1)) / 19
    ctrl = {'symbol': symbol}
    obs = []
    for sp in species:
        s = seqs.get(sp)
        if s is None:
            continue
        for i, c in enumerate(s):
            if c in IDX and c != anc[i] and plddt[i] >= 70:
                obs.append((dd_all[i, IDX[c]], alt_mean[i]))
    if obs:
        o = np.array(obs)
        ctrl.update({'n_substitution_events': len(o), 'mean_observed_ddG': o[:, 0].mean(),
                     'mean_alternative_ddG_same_sites': o[:, 1].mean()})
    el = plddt >= 70
    ctrl.update({'alt_abs_ddG_buried': float(np.abs(dd_all[el & (rsa < 0.25)]).sum() / max(1, 19 * (el & (rsa < 0.25)).sum())),
                 'alt_abs_ddG_exposed': float(np.abs(dd_all[el & (rsa >= 0.25)]).sum() / max(1, 19 * (el & (rsa >= 0.25)).sum()))})
    for cls in ['helix_region', 'strand_region']:
        m = el & (ss == cls) & (a_idx != IDX['P']) & (a_idx != IDX['G'])
        if m.any():
            ctrl[f'ddG_to_P_{cls}'] = float(dd_all[m, IDX['P']].mean())
            ctrl[f'ddG_to_G_{cls}'] = float(dd_all[m, IDX['G']].mean())
            ctrl[f'ddG_mean_alt_{cls}'] = float(alt_mean[m].mean())
    return rows, ctrl


def main():
    sp = pd.read_csv(DATA / 'species_final.csv')
    species = sorted(sp[sp.eligible & sp.placental].ensembl_name)
    genes = pd.read_csv(RESULTS / 'thermompnn_summary.csv')
    genes = genes[genes.status == 'ok'].symbol.tolist()
    with ProcessPoolExecutor(max_workers=16) as pool:
        res = list(pool.map(gene_metrics, [(g, species) for g in genes], chunksize=8))
    metrics = pd.DataFrame([r for rows, _ in res for r in rows])
    controls = pd.DataFrame([c for _, c in res])
    metrics.to_csv(RESULTS / 'species_gene_metrics.csv.gz', index=False)
    controls.to_csv(RESULTS / 'predictor_controls_by_gene.csv', index=False)
    c = controls
    summary = {
        'genes': len(c),
        'C1_mean_observed_ddG': float(np.average(c.mean_observed_ddG.dropna(), weights=c.n_substitution_events.dropna())),
        'C1_mean_alternative_ddG': float(np.average(c.mean_alternative_ddG_same_sites.dropna(), weights=c.n_substitution_events.dropna())),
        'C1_fraction_genes_observed_below_alternatives': float((c.mean_observed_ddG < c.mean_alternative_ddG_same_sites).mean()),
        'C2_mean_abs_alt_buried': float(c.alt_abs_ddG_buried.mean()), 'C2_mean_abs_alt_exposed': float(c.alt_abs_ddG_exposed.mean()),
        'C2_fraction_genes_buried_gt_exposed': float((c.alt_abs_ddG_buried > c.alt_abs_ddG_exposed).mean()),
    }
    for cls in ['helix_region', 'strand_region']:
        summary[f'C3_{cls}_to_P'] = float(c[f'ddG_to_P_{cls}'].mean())
        summary[f'C3_{cls}_to_G'] = float(c[f'ddG_to_G_{cls}'].mean())
        summary[f'C3_{cls}_mean_alt'] = float(c[f'ddG_mean_alt_{cls}'].mean())
    pd.Series(summary).to_json(RESULTS / 'predictor_controls_summary.json', indent=2)
    print(pd.Series(summary).to_string())
    print('species-gene pairs', len(metrics), '| species', metrics.species.nunique(), '| genes', metrics.symbol.nunique())


if __name__ == '__main__':
    main()
