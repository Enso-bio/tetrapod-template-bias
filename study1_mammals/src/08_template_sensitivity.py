"""Template sensitivity (plan section 7): mouse AlphaFold backbone instead of human for 200 random eligible genes.

The ancestral sequence is threaded onto the mouse backbone at positions aligned to human; species metrics are then
recomputed on the sites scored with both templates and compared with the human-template values.
"""
import io
from urllib.parse import urlencode

import numpy as np
import pandas as pd
from Bio import Phylo

from common import DATA, RAW, RESULTS, SEED, fetch, fetch_json
from thermo import AA20, load_model, parse, site_energies
from pgls import PGLS, summarize, tree_covariance

IDX = {a: i for i, a in enumerate(AA20)}
N_GENES = 200


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


def mouse_index_map(human, mouse_seq):
    """Human position -> mouse residue index, with the same semi-global aligner used for the ortholog alignments."""
    from Bio.Align import PairwiseAligner, substitution_matrices
    al = PairwiseAligner(mode='global', substitution_matrix=substitution_matrices.load('BLOSUM62'),
                         open_gap_score=-10, extend_gap_score=-0.5, end_insertion_score=0, end_deletion_score=0)
    clean = lambda s: ''.join(c if c in AA20 + 'X' else 'X' for c in s.replace('U', 'C'))
    a = al.align(clean(human), clean(mouse_seq))[0]
    m = {}
    for (h0, h1), (o0, o1) in zip(*a.aligned):
        for i, j in zip(range(h0, h1), range(o0, o1)):
            m[i] = j
    return m


def main():
    rng = np.random.default_rng(SEED + 8)
    summary = pd.read_csv(RESULTS / 'thermompnn_summary.csv')
    genes_ok = summary[summary.status == 'ok'].symbol.tolist()
    sel = pd.read_csv(DATA / 'gene_selection_audit.csv').set_index('symbol')
    o2o = pd.read_csv(DATA / 'one2one_orthologs.csv.gz')
    o2o = o2o[o2o.species == 'mus_musculus'].set_index('human_gene')
    q = {'query': '(organism_id:10090) AND (reviewed:true)', 'format': 'tsv', 'fields': 'accession,sequence'}
    mouse_uni = pd.read_csv(io.BytesIO(fetch('https://rest.uniprot.org/uniprotkb/stream?' + urlencode(q),
                                             RAW / 'uniprot_mouse_swissprot.tsv', timeout=600)), sep='\t')
    seq_to_acc = mouse_uni.drop_duplicates('Sequence', keep=False).set_index('Sequence').Entry
    import gzip
    mpep = {}
    with gzip.open(RAW / 'pep' / 'mus_musculus.pep.all.fa.gz', 'rt') as f:
        pid, buf = None, []
        for line in f:
            if line.startswith('>'):
                if pid:
                    mpep[pid] = ''.join(buf)
                pid, buf = line[1:].split()[0].split('.')[0], []
            else:
                buf.append(line.strip())
        mpep[pid] = ''.join(buf)

    order = rng.permutation(genes_ok)
    model = load_model()
    sp = pd.read_csv(DATA / 'species_final.csv')
    species = sorted(sp[sp.eligible & sp.placental].ensembl_name)
    rows, used = [], []
    for g in order:
        if len(used) >= N_GENES:
            break
        hg = sel.loc[g, 'human_gene']
        if hg not in o2o.index:
            continue
        mseq = mpep.get(o2o.loc[hg, 'protein'])
        acc = seq_to_acc.get(mseq) if mseq else None
        if acc is None:
            continue
        info = fetch_json(f'https://alphafold.ebi.ac.uk/api/prediction/{acc}', RAW / 'alphafold_api' / f'{acc}.json', ok_404=True)
        full = [m for m in (info or []) if m.get('sequence', m.get('uniprotSequence')) == mseq]
        if not full:
            continue
        pdb = fetch(max(full, key=lambda m: m.get('latestVersion', 0))['pdbUrl'], DATA / 'alphafold_mouse' / f'{acc}.pdb')
        z = np.load(DATA / 'energies' / f'{g}.npz')
        anc, human, E_h, P, plddt_h = str(z['ancestor']), str(z['human']), z['e_anc'], z['probs'], z['plddt']
        d = parse(DATA / 'alphafold_mouse' / f'{acc}.pdb')
        if len(d['seq_chain_A']) != len(mseq):
            continue
        hmap = mouse_index_map(human, mseq)
        threaded = list(d['seq_chain_A'])
        for i, j in hmap.items():
            threaded[j] = anc[i]
        E_m_full = site_energies(model, d, ''.join(threaded))
        sites = np.array(sorted(i for i in hmap if plddt_h[i] >= 70))
        E_m = np.full_like(E_h, np.nan)
        E_m[sites] = E_m_full[[hmap[i] for i in sites]]
        seqs = read_fasta(DATA / 'alignments' / f'{g}.fasta')
        for s in species:
            sq = seqs.get(s)
            if sq is None:
                continue
            for i in sites:
                c = sq[i].replace('U', 'C')
                if c not in IDX:
                    continue
                for tag, E in (('human', E_h), ('mouse', E_m)):
                    ddg = E[i, IDX[c]] - (P[i] * E[i]).sum()
                    rows.append((g, s, tag, ddg, 1 - P[i, IDX[c]]))
        used.append(g)
    df = pd.DataFrame(rows, columns=['symbol', 'species', 'template', 'ddg', 'k'])
    agg = df.groupby(['species', 'template']).apply(lambda x: x.ddg.sum() / x.k.sum(), include_groups=False).unstack()
    site_corr = float(np.corrcoef(df[df.template == 'human'].ddg, df[df.template == 'mouse'].ddg)[0, 1])
    sp2 = sp.set_index('ensembl_name').loc[agg.index].reset_index()
    tree = Phylo.read(RAW / 'mcc_node.tre', 'nexus')
    model_p = PGLS(tree_covariance(tree, sp2.tree_label.tolist()))
    y = np.log10(sp2.max_lifespan_yr.to_numpy())
    res = {'genes_used': len(used), 'species': len(agg), 'site_ddG_correlation_human_vs_mouse_template': site_corr,
           'species_M_correlation': float(np.corrcoef(agg.human, agg.mouse)[0, 1])}
    for tag in ('human', 'mouse'):
        m = agg[tag].to_numpy()
        X = np.column_stack([np.ones(len(y)), np.log10(sp2.adult_mass_g), (m - m.mean()) / m.std(ddof=1)])
        res[f'P1_like_{tag}_template'] = summarize(model_p.fit(y, X), 2)
    pd.Series(used).to_csv(RESULTS / 'template_sensitivity_genes.csv', index=False, header=['symbol'])
    agg.to_csv(RESULTS / 'template_sensitivity_species_M.csv')
    pd.Series(res).to_json(RESULTS / 'template_sensitivity.json', indent=2, default_handler=str)
    print(res)


if __name__ == '__main__':
    main()
