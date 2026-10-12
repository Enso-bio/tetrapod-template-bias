"""Template sensitivity (plan section 6): coelacanth backbone (all species) and human backbone (reptiles only),
for 300 random eligible genes, on sites shared with the zebrafish template."""
import io
import json
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from Bio import Phylo
from Bio.Align import PairwiseAligner, substitution_matrices

from common import DATA, RAW, RESULTS, SEED, fetch, fetch_json
from pgls import PGLS, summarize, tree_covariance
from thermo import AA20, load_model, parse, site_energies

IDX = {a: i for i, a in enumerate(AA20)}
TEMPLATES = {'latimeria_chalumnae': 7897, 'homo_sapiens': 9606}
N_GENES, WINDOW, MAX_DIFF = 300, 10, 5
REPTILES = {'Testudines', 'Crocodylia', 'Lepidosauria'}


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


def index_map(a, b):
    al = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                         open_gap_score=-10, extend_gap_score=-0.5, end_insertion_score=0, end_deletion_score=0)
    clean = lambda s: ''.join(c if c in AA20 else ('C' if c == 'U' else 'X') for c in s)
    x = al.align(clean(a), clean(b))[0]
    return {i: j for (a0, a1), (b0, b1) in zip(*x.aligned) for i, j in zip(range(a0, a1), range(b0, b1))}


def window_mask(diff):
    mask = np.zeros(len(diff), bool)
    c = np.concatenate([[0], np.cumsum(diff)])
    for w in range(0, max(1, len(diff) - WINDOW + 1)):
        if c[min(w + WINDOW, len(diff))] - c[w] >= MAX_DIFF:
            mask[w:w + WINDOW] = True
    return mask


def main():
    rng = np.random.default_rng(SEED + 8)
    genes = pd.read_csv(RESULTS / 'thermompnn_summary.csv').query("status == 'ok'").human_gene.tolist()
    order = list(rng.permutation(genes))
    seq_to_acc = {}
    for name, tax in TEMPLATES.items():
        q = {'query': f'(organism_id:{tax})', 'format': 'tsv', 'fields': 'accession,reviewed,sequence'}
        up = pd.read_csv(io.BytesIO(fetch('https://rest.uniprot.org/uniprotkb/stream?' + urlencode(q),
                                          RAW / f'uniprot_{name}.tsv', timeout=900)), sep='\t')
        up = up.sort_values(['Reviewed', 'Entry'], ascending=[False, True]).drop_duplicates('Sequence')
        seq_to_acc[name] = up.set_index('Sequence').Entry
    sp = pd.read_csv(DATA / 'species_final.csv')
    sp = sp[sp.eligible].sort_values('binomial').reset_index(drop=True)
    species = sp.binomial.tolist()
    model = load_model()
    rows, used = [], []
    for g in order:
        if len(used) >= N_GENES:
            break
        raw = read_fasta(DATA / 'sequences' / f'{g}.fasta')
        z = np.load(DATA / 'energies' / f'{g}.npz')
        anc, zf, E_zf, P, conf = str(z['ancestor']), str(z['template']), z['e_anc'], z['probs'], z['plddt']
        E_alt, ok_alt = {}, {}
        for name in TEMPLATES:
            tseq = raw.get(name)
            acc = seq_to_acc[name].get(tseq) if tseq else None
            if acc is None:
                break
            info = fetch_json(f'https://alphafold.ebi.ac.uk/api/prediction/{acc}', RAW / 'alphafold_api' / f'{acc}.json', ok_404=True)
            full = [m for m in (info or []) if m.get('sequence', m.get('uniprotSequence')) == tseq]
            if not full:
                break
            pdb = fetch(max(full, key=lambda m: m.get('latestVersion', 0))['pdbUrl'], DATA / f'alphafold_{name}' / f'{acc}.pdb')
            d = parse(DATA / f'alphafold_{name}' / f'{acc}.pdb')
            if len(d['seq_chain_A']) != len(tseq):
                break
            m = index_map(zf, tseq)
            threaded = list(d['seq_chain_A'])
            for i, j in m.items():
                threaded[j] = anc[i]
            e = site_energies(model, d, ''.join(threaded))
            t_conf = np.array([float(l[60:66]) for l in open(DATA / f'alphafold_{name}' / f'{acc}.pdb', encoding='utf-8')
                               if l.startswith('ATOM') and l[12:16].strip() == 'CA'])
            E = np.full((len(zf), 20), np.nan)
            ok = np.zeros(len(zf), bool)
            for i, j in m.items():
                E[i], ok[i] = e[j], t_conf[j] >= 70
            E_alt[name], ok_alt[name] = E, ok
        if len(E_alt) < len(TEMPLATES):
            continue
        sites = (conf >= 70) & ok_alt['latimeria_chalumnae'] & ok_alt['homo_sapiens']
        aln = read_fasta(DATA / 'alignments' / f'{g}.fasta')
        exp = {k: np.nansum(P * np.nan_to_num(v), axis=1) for k, v in
               {'zebrafish': E_zf, **E_alt}.items()}
        Es = {'zebrafish': E_zf, **E_alt}
        for s in species:
            q = aln.get(s)
            if q is None:
                continue
            code = np.array([IDX.get(c.replace('U', 'C'), -1) for c in q])
            present = code >= 0
            masked = window_mask((present & (np.array(list(q)) != np.array(list(anc)))).astype(int))
            msk = present & sites & ~masked
            ci = np.where(present, code, 0)
            r = np.arange(len(q))
            rec = {'human_gene': g, 'species': s, 'K': float((1 - P[r, ci])[msk].sum())}
            for k, E in Es.items():
                rec[f'D_{k}'] = float((E[r, ci] - exp[k])[msk].sum())
            rows.append(rec)
        used.append(g)
    df = pd.DataFrame(rows)
    agg = df.groupby('species')[['K', 'D_zebrafish', 'D_latimeria_chalumnae', 'D_homo_sapiens']].sum()
    M = pd.DataFrame({k: agg[f'D_{k}'] / agg.K for k in ['zebrafish', 'latimeria_chalumnae', 'homo_sapiens']}).loc[species]
    M.to_csv(RESULTS / 'template_sensitivity_species_M.csv')
    tree = Phylo.read(DATA / 'timetree_upgma.nwk', 'newick')
    zf_ = lambda x: (x - x.mean()) / x.std(ddof=1)
    res = {'genes_used': len(used), 'species': len(species),
           'r_species_zebrafish_coelacanth': float(np.corrcoef(M.zebrafish, M.latimeria_chalumnae)[0, 1])}
    for label, keep, tpls in [('all_species', np.ones(len(sp), bool), ['zebrafish', 'latimeria_chalumnae']),
                              ('reptiles_only', sp.group.isin(REPTILES).to_numpy(), ['zebrafish', 'homo_sapiens'])]:
        sub = sp[keep].reset_index(drop=True)
        model_p = PGLS(tree_covariance(tree, sub.binomial.tolist()))
        y = np.log10(sub.max_lifespan_yr.to_numpy())
        for t in tpls:
            m = M[t].to_numpy()[keep]
            X = np.column_stack([np.ones(len(y)), np.log10(sub.mass_g), zf_(m)])
            res[f'P1_like_{label}_{t}'] = summarize(model_p.fit(y, X), 2)
        if label == 'reptiles_only':
            res['r_species_reptiles_zebrafish_human'] = float(np.corrcoef(M.zebrafish[keep], M.homo_sapiens[keep])[0, 1])
    turt = (sp.group == 'Testudines').to_numpy()
    for t in M.columns:
        v = M[t].to_numpy()
        res[f'turtle_offset_SD_{t}'] = float((v[turt].mean() - v[~turt].mean()) / v.std(ddof=1))
    pd.Series(used).to_csv(RESULTS / 'template_sensitivity_genes.csv', index=False, header=['human_gene'])
    (RESULTS / 'template_sensitivity.json').write_text(json.dumps(res, indent=2, default=float))
    print(json.dumps(res, indent=1, default=float))


if __name__ == '__main__':
    main()
