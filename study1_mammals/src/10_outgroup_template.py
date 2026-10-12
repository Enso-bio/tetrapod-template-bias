"""Stability metrics on an outgroup backbone equidistant from all placental species.

Usage: python 10_outgroup_template.py [ensembl_name taxon_id role]
  default: monodelphis_domestica 13616 post_hoc_E4   (opossum; plan amendment E4, post hoc)
  E5:      ornithorhynchus_anatinus 9258 E5_template (platypus; energies for the confirmatory tests in 12_*.py)

Template: AlphaFold DB model of the outgroup one-to-one ortholog with a sequence identical to a UniProt entry. The
reconstructed placental ancestor is threaded at positions aligned to human; sites need pLDDT >= 70 in the human and
in the template model. The same window mask, metrics and PGLS procedures as the main analysis are applied.
"""
import gzip
import io
import json
import sys
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Align import PairwiseAligner, substitution_matrices

from common import DATA, RAW, RESULTS, SEED, fetch, fetch_json
from pgls import PGLS, bootstrap_lrt, summarize, tree_covariance
from thermo import AA20, load_model, parse, site_energies

TEMPLATE, TAXON, ROLE = (sys.argv[1], int(sys.argv[2]), sys.argv[3]) if len(sys.argv) > 3 else \
    ('monodelphis_domestica', 13616, 'post_hoc_E4')
IDX = {a: i for i, a in enumerate(AA20)}
WINDOW, MAX_DIFF = 10, 5
OUT = DATA / f'energies_{TEMPLATE}'
OUT.mkdir(exist_ok=True)


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


def read_pep(path):
    out, pid, buf = {}, None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if pid:
                    out[pid] = ''.join(buf)
                pid, buf = line[1:].split()[0].split('.')[0], []
            else:
                buf.append(line.strip())
    out[pid] = ''.join(buf)
    return out


def index_map(human, other):
    al = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                         open_gap_score=-10, extend_gap_score=-0.5, end_insertion_score=0, end_deletion_score=0)
    clean = lambda s: ''.join(c if c in AA20 else ('C' if c == 'U' else 'X') for c in s)
    a = al.align(clean(human), clean(other))[0]
    return {i: j for (h0, h1), (o0, o1) in zip(*a.aligned) for i, j in zip(range(h0, h1), range(o0, o1))}


def window_mask(diff):
    mask = np.zeros(len(diff), bool)
    c = np.concatenate([[0], np.cumsum(diff)])
    for w in range(0, max(1, len(diff) - WINDOW + 1)):
        if c[min(w + WINDOW, len(diff))] - c[w] >= MAX_DIFF:
            mask[w:w + WINDOW] = True
    return mask


def energies():
    genes = pd.read_csv(RESULTS / 'thermompnn_summary.csv').query("status == 'ok'").symbol.tolist()
    sel = pd.read_csv(DATA / 'gene_selection_audit.csv').set_index('symbol')
    o2o = pd.read_csv(DATA / 'one2one_orthologs.csv.gz')
    o2o = o2o[o2o.species == TEMPLATE].set_index('human_gene')
    q = {'query': f'(organism_id:{TAXON})', 'format': 'tsv', 'fields': 'accession,reviewed,sequence'}
    uni = pd.read_csv(io.BytesIO(fetch('https://rest.uniprot.org/uniprotkb/stream?' + urlencode(q),
                                       RAW / f'uniprot_{TEMPLATE}.tsv', timeout=900)), sep='\t')
    uni = uni.sort_values(['Reviewed', 'Entry'], ascending=[False, True]).drop_duplicates('Sequence')
    seq_to_acc = uni.set_index('Sequence').Entry
    pep = read_pep(RAW / 'pep' / f'{TEMPLATE}.pep.all.fa.gz')
    model = load_model()
    audit = []
    for g in genes:
        path = OUT / f'{g}.npz'
        if path.exists():
            audit.append({'symbol': g, 'status': 'ok'})
            continue
        hg = sel.loc[g, 'human_gene']
        tseq = pep.get(o2o.loc[hg, 'protein']) if hg in o2o.index else None
        acc = seq_to_acc.get(tseq) if tseq else None
        if acc is None:
            audit.append({'symbol': g, 'status': 'no_identical_uniprot_entry'})
            continue
        info = fetch_json(f'https://alphafold.ebi.ac.uk/api/prediction/{acc}', RAW / 'alphafold_api' / f'{acc}.json', ok_404=True)
        full = [m for m in (info or []) if m.get('sequence', m.get('uniprotSequence')) == tseq]
        if not full:
            audit.append({'symbol': g, 'status': 'no_alphafold_model', 'accession': acc})
            continue
        pdb = DATA / f'alphafold_{TEMPLATE}' / f'{acc}.pdb'
        fetch(max(full, key=lambda m: m.get('latestVersion', 0))['pdbUrl'], pdb)
        d = parse(pdb)
        if len(d['seq_chain_A']) != len(tseq):
            audit.append({'symbol': g, 'status': 'pdb_length_mismatch', 'accession': acc})
            continue
        z = np.load(DATA / 'energies' / f'{g}.npz')
        anc, human = str(z['ancestor']), str(z['human'])
        m = index_map(human, tseq)
        threaded = list(d['seq_chain_A'])
        for i, j in m.items():
            threaded[j] = anc[i]
        e_full = site_energies(model, d, ''.join(threaded))
        t_plddt = np.array([float(l[60:66]) for l in open(pdb, encoding='utf-8')
                            if l.startswith('ATOM') and l[12:16].strip() == 'CA'])
        E = np.full((len(human), 20), np.nan)
        ok = np.zeros(len(human), bool)
        for i, j in m.items():
            E[i] = e_full[j]
            ok[i] = t_plddt[j] >= 70
        np.savez_compressed(path, e=E, template_ok=ok, accession=acc)
        audit.append({'symbol': g, 'status': 'ok', 'accession': acc, 'mapped_sites': len(m)})
    pd.DataFrame(audit).to_csv(RESULTS / f'template_{TEMPLATE}_audit.csv', index=False)
    return pd.DataFrame(audit)


def metrics(species):
    rows = []
    for f in sorted(OUT.glob('*.npz')):
        g = f.stem
        t = np.load(f)
        z = np.load(DATA / 'energies' / f'{g}.npz')
        E, P, anc, plddt = t['e'], z['probs'], str(z['ancestor']), z['plddt']
        base_ok = t['template_ok'] & (plddt >= 70)
        expected = np.nansum(P * E, axis=1)
        seqs = read_fasta(DATA / 'alignments' / f'{g}.fasta')
        for sp in species:
            s = seqs.get(sp)
            if s is None:
                continue
            code = np.array([IDX.get(c.replace('U', 'C'), -1) for c in s])
            present = code >= 0
            masked = window_mask((present & (np.array(list(s)) != np.array(list(anc)))).astype(int))
            m = present & base_ok & ~masked
            ci = np.where(present, code, 0)
            ddg = E[np.arange(len(s)), ci] - expected
            k = 1 - P[np.arange(len(s)), ci]
            rows.append({'symbol': g, 'species': sp, 'D': float(ddg[m].sum()), 'K': float(k[m].sum()), 'L': int(m.sum())})
    return pd.DataFrame(rows)


def main():
    audit = energies()
    print(audit.status.value_counts().to_dict(), flush=True)
    sp = pd.read_csv(DATA / 'species_final.csv')
    sp = sp[sp.eligible & sp.placental].sort_values('ensembl_name').reset_index(drop=True)
    met = metrics(sp.ensembl_name.tolist())
    met.to_csv(RESULTS / f'species_gene_metrics_{TEMPLATE}.csv.gz', index=False)
    genes = pd.read_csv(RESULTS / 'asr_summary.csv').set_index('symbol')
    symbols = sorted(met.symbol.unique())
    D = met.pivot(index='species', columns='symbol', values='D').reindex(index=sp.ensembl_name, columns=symbols).to_numpy()
    K = met.pivot(index='species', columns='symbol', values='K').reindex(index=sp.ensembl_name, columns=symbols).to_numpy()
    gset = genes.gene_set.reindex(symbols).to_numpy()
    agg = lambda cols: np.nansum(D[:, cols], 1) / np.nansum(K[:, cols], 1)
    allg = np.ones(len(symbols), bool)
    gen, bg = gset == 'genage', gset == 'background'
    tree = Phylo.read(RAW / 'mcc_node.tre', 'nexus')
    y = np.log10(sp.max_lifespan_yr.to_numpy())
    lm = np.log10(sp.adult_mass_g.to_numpy())
    zf = lambda x: (x - x.mean()) / x.std(ddof=1)
    model = PGLS(tree_covariance(tree, sp.tree_label.tolist()))
    out = {'template': TEMPLATE, 'role': ROLE, 'post_hoc': ROLE.startswith('post_hoc'), 'n_species': len(sp), 'n_genes': int(allg.sum()), 'n_genage': int(gen.sum())}
    M_all, M_gen = agg(allg), agg(gen)
    for name, m, seed in [('P1_posthoc', M_all, SEED + 41), ('P2_selfcontained_posthoc', M_gen, SEED + 42)]:
        X = np.column_stack([np.ones(len(y)), lm, zf(m)])
        fit, b = bootstrap_lrt(model, y, X, 2, 999, seed)
        out[name] = {**summarize(fit, 2), **b}
    rng = np.random.default_rng(SEED + 43)
    bg_idx = np.where(bg)[0]
    null = []
    for _ in range(1000):
        cols = np.zeros(len(symbols), bool)
        cols[rng.choice(bg_idx, gen.sum(), replace=False)] = True
        null.append(model.fit(y, np.column_stack([np.ones(len(y)), lm, zf(agg(cols))]))['beta'][2])
    null = np.array(null)
    bo = out['P2_selfcontained_posthoc']['beta']
    out['P2_competitive_posthoc'] = {'beta_genage': bo, 'null_mean': float(null.mean()), 'null_sd': float(null.std()),
                                     'p_empirical_two_sided': float((np.sum(np.abs(null) >= abs(bo)) + 1) / 1001)}
    for label, keep in [('exclude_primates', sp.order != 'Primates'), ('exclude_human', sp.ensembl_name != 'homo_sapiens')]:
        k = keep.to_numpy()
        mk = PGLS(tree_covariance(tree, sp.tree_label[k].tolist()))
        out[f'P1_posthoc_{label}'] = summarize(mk.fit(y[k], np.column_stack([np.ones(k.sum()), lm[k], zf(M_all[k])])), 2)
    main_st = pd.read_csv(RESULTS / 'species_metrics.csv').set_index('ensembl_name').loc[sp.ensembl_name]
    st = pd.DataFrame({'ensembl_name': sp.ensembl_name, 'order': sp.order, 'M_all_outgroup': M_all, 'M_genage_outgroup': M_gen,
                       'M_all_human_template': main_st.M_all.to_numpy(), 'LQ': sp.LQ})
    st.to_csv(RESULTS / f'species_metrics_{TEMPLATE}.csv', index=False)
    zs = zf(M_all)
    out['by_order_mean_z'] = st.assign(z=zs).groupby('order').z.mean().round(3).to_dict()
    out['human_rank_of_57'] = int(pd.Series(M_all).rank().iloc[list(sp.ensembl_name).index('homo_sapiens')])
    out['lowest_5_species'] = st.sort_values('M_all_outgroup').ensembl_name.head(5).tolist()
    out['highest_5_species'] = st.sort_values('M_all_outgroup').ensembl_name.tail(5).tolist()
    out['corr_with_human_template'] = float(np.corrcoef(M_all, st.M_all_human_template)[0, 1])
    out['nearly_neutral_M_on_log10mass_posthoc'] = summarize(
        model.fit(zf(M_all), np.column_stack([np.ones(len(y)), lm])), 1)
    name = f'posthoc_{TEMPLATE}_template.json' if ROLE.startswith('post_hoc') else f'template_{TEMPLATE}_summary.json'
    (RESULTS / name).write_text(json.dumps(out, indent=2, default=float))
    np.save(RESULTS / f'P2_null_betas_{TEMPLATE}.npy', null)
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
