"""POST HOC amendment E6: human-independent alignment and sites with the platypus template.

Genes: the 1,341 common genes of E5. For each gene the Compara 1:1 proteins of all analysed species (human included as
one more species) are aligned with MAFFT and the columns with a platypus residue are kept. Eligible sites: pLDDT >= 70
in the platypus model. The placental ancestor is reconstructed on this alignment and threaded onto the platypus model.
No step uses the human sequence or model as a reference. Variant 'human_sites' restricts the same metric to columns
whose aligned human residue has pLDDT >= 70 in the human model (the site definition of E4/E5).
"""
import gzip
import importlib
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
from Bio import Phylo

from common import DATA, RAW, RESULTS, TOOLS
from pgls import PGLS, bootstrap_lrt, summarize, tree_covariance
from thermo import load_model, parse, site_energies

asr = importlib.import_module('04_asr')
e5 = importlib.import_module('12_confirm_template')

MAFFT = TOOLS / 'mafft-win' / 'mafft.bat'
IQTREE = asr.IQTREE
TEMPLATE = 'ornithorhynchus_anatinus'
WORK = DATA / 'msa_platypus'
SEQ, ALN, ASR, ENERGY = WORK / 'sequences', WORK / 'alignments', WORK / 'asr', WORK / 'energies'
for folder in (SEQ, ALN, ASR, ENERGY):
    folder.mkdir(parents=True, exist_ok=True)
IQ_AA = 'ARNDCQEGHILKMFPSTWYV'
AA20 = 'ACDEFGHIKLMNPQRSTVWY'
TO_AA20 = [IQ_AA.index(a) for a in AA20]
IDX = {a: i for i, a in enumerate(AA20)}
SEED_BOOT, SEED_P1, B_GENES = 20261070, 20261071, 1000


def clean(s):
    s = s.rstrip('*').upper().replace('U', 'C')
    return ''.join(c if c in AA20 else 'X' for c in s)


def read_pep(path, wanted):
    out, pid, buf = {}, None, []
    with gzip.open(path, 'rt') as f:
        for line in f:
            if line.startswith('>'):
                if pid in wanted:
                    out[pid] = ''.join(buf)
                pid, buf = line[1:].split()[0].split('.')[0], []
            else:
                buf.append(line.strip())
    if pid in wanted:
        out[pid] = ''.join(buf)
    return out


def plddt_ca(pdb):
    return np.array([float(l[60:66]) for l in open(pdb, encoding='utf-8')
                     if l.startswith('ATOM') and l[12:16].strip() == 'CA'])


def write_sequences(genes, species, sel):
    """One FASTA per gene with the Compara 1:1 protein of every species that has one."""
    o2o = pd.read_csv(DATA / 'one2one_orthologs.csv.gz')
    hg = sel.set_index('symbol').human_gene.reindex(genes)
    o2o = o2o[o2o.human_gene.isin(set(hg)) & o2o.species.isin(species)]
    human = sel.set_index('symbol').human_protein.reindex(genes)
    wanted = {s: set(o2o[o2o.species == s].protein) for s in species}
    wanted['homo_sapiens'] = set(human)
    seqs = {}
    for s in species:
        seqs[s] = read_pep(RAW / 'pep' / f'{s}.pep.all.fa.gz', wanted[s])
    prot = o2o.set_index(['human_gene', 'species']).protein
    for g in genes:
        path = SEQ / f'{g}.fasta'
        if path.exists():
            continue
        with open(path, 'w', encoding='utf-8') as f:
            for s in species:
                pid = human[g] if s == 'homo_sapiens' else prot.get((hg[g], s))
                if pid is not None and pid in seqs[s]:
                    f.write(f'>{s}\n{clean(seqs[s][pid])}\n')


def align(g):
    """MAFFT, then keep the platypus columns; also record the human residue index of each kept column."""
    path = ALN / f'{g}.fasta'
    if path.exists():
        return g, 'ok'
    p = subprocess.run([str(MAFFT), '--quiet', '--auto', '--thread', '1', str(SEQ / f'{g}.fasta')],
                       capture_output=True, text=True)
    msa = parse_fasta(p.stdout)
    if TEMPLATE not in msa or 'homo_sapiens' not in msa:
        return g, 'mafft_failed'
    cols = [i for i, c in enumerate(msa[TEMPLATE]) if c != '-']
    hpos = np.cumsum([c != '-' for c in msa['homo_sapiens']]) - 1
    human_idx = np.array([hpos[i] if msa['homo_sapiens'][i] != '-' else -1 for i in cols])
    np.save(ALN / f'{g}.human_idx.npy', human_idx)
    with open(path, 'w', encoding='utf-8') as f:
        for s, q in msa.items():
            f.write(f'>{s}\n{"".join(q[i] for i in cols).upper()}\n')
    return g, 'ok'


def parse_fasta(text):
    seqs, name = {}, None
    for line in text.splitlines():
        line = line.strip()
        if line.startswith('>'):
            name = line[1:]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: ''.join(v) for k, v in seqs.items()}


def run_asr(args):
    g, names, placentals, nwk = args
    folder = ASR / g
    folder.mkdir(exist_ok=True)
    state = folder / f'{g}.state'
    if not state.exists():
        seqs = asr.read_fasta(ALN / f'{g}.fasta')
        with open(folder / 'aln.fasta', 'w', encoding='utf-8') as f:
            for s in names:
                f.write(f'>{s}\n{seqs[s]}\n')
        (folder / 'tree.nwk').write_text(nwk)
        p = subprocess.run([str(IQTREE), '-s', 'aln.fasta', '-te', 'tree.nwk', '-m', 'LG+G4', '-asr', '-nt', '1',
                            '-pre', g, '-quiet', '-redo', '-seed', '20261008'], cwd=folder, capture_output=True, text=True)
        if p.returncode != 0 or not state.exists():
            return {'symbol': g, 'status': 'iqtree_failed', 'stderr': (p.stderr or p.stdout)[-300:]}
    node = asr.placental_ancestor_node(folder / f'{g}.treefile', placentals)
    st = pd.read_csv(state, sep='\t', comment='#')
    probs = st[st.Node == node].sort_values('Site')[[f'p_{a}' for a in IQ_AA]].to_numpy()[:, TO_AA20]
    return {'symbol': g, 'status': 'ok', 'node': node, 'probs': probs, 'mean_max_posterior': float(probs.max(1).mean())}


def gene_metrics(g, species, probs, pdb, human_plddt):
    """D, K, L per species for the main metric and the human_sites variant."""
    z = np.load(ENERGY / f'{g}.npz')
    E, anc = z['e_anc'], str(z['ancestor'])
    ok = plddt_ca(pdb) >= 70
    hidx = np.load(ALN / f'{g}.human_idx.npy')
    human_ok = np.array([i >= 0 and human_plddt[i] >= 70 for i in hidx])
    expected = (probs * E).sum(axis=1)
    seqs = asr.read_fasta(ALN / f'{g}.fasta')
    rows = []
    for s in species:
        q = seqs.get(s)
        if q is None:
            continue
        code = np.array([IDX.get(c, -1) for c in q])
        present = code >= 0
        masked = e5.window_mask((present & (np.array(list(q)) != np.array(list(anc)))).astype(int))
        ci = np.where(present, code, 0)
        r = np.arange(len(q))
        ddg, k = E[r, ci] - expected, 1 - probs[r, ci]
        rec = {'symbol': g, 'species': s}
        for name, m in [('main', present & ok & ~masked), ('human_sites', present & ok & human_ok & ~masked)]:
            rec[f'D_{name}'], rec[f'K_{name}'], rec[f'L_{name}'] = float(ddg[m].sum()), float(k[m].sum()), int(m.sum())
        rows.append(rec)
    return rows


def offset(M, prim):
    return (M[prim].mean() - M[~prim].mean()) / M.std(ddof=1)


def main():
    sp_all = pd.read_csv(DATA / 'species_final.csv')
    sp_all = sp_all[sp_all.eligible]
    species = sorted(sp_all.ensembl_name)
    placentals = sorted(sp_all[sp_all.placental].ensembl_name)
    outgroups = sorted(set(species) - set(placentals))
    sets = [set(f.stem for f in (DATA / f'energies_{t}').glob('*.npz')) for t in e5.TEMPLATES.values()]
    genes = sorted(set.intersection(*sets))
    sel = pd.read_csv(DATA / 'gene_selection_audit.csv')
    sel = sel[sel.symbol.isin(genes)]
    write_sequences(genes, species, sel)
    with ThreadPoolExecutor(max_workers=14) as pool:
        st = dict(pool.map(align, genes))
    print('MAFFT', pd.Series(st).value_counts().to_dict(), flush=True)

    # Coverage on platypus-eligible sites, gene filter (>= 60% of placentals), per-gene species lists
    acc = {g: str(np.load(DATA / f'energies_{TEMPLATE}' / f'{g}.npz')['accession']) for g in genes}
    pdbs = {g: DATA / f'alphafold_{TEMPLATE}' / f'{acc[g]}.pdb' for g in genes}
    audit, tasks = [], []
    labels = dict(zip(sp_all.tree_label, sp_all.ensembl_name))
    base_tree = asr.prune_tree(RAW / 'mcc_node.tre', labels)
    for g in genes:
        if st[g] != 'ok':
            audit.append({'symbol': g, 'status': st[g]})
            continue
        seqs = asr.read_fasta(ALN / f'{g}.fasta')
        model_seq = parse(pdbs[g])['seq_chain_A'].replace('U', 'C')
        if len(model_seq) != len(seqs[TEMPLATE]) or any(a != b and b != 'X' for a, b in zip(model_seq, seqs[TEMPLATE])):
            audit.append({'symbol': g, 'status': 'template_mismatch'})
            continue
        el = plddt_ca(pdbs[g]) >= 70
        cov = {s: float(np.mean([c not in '-X' for c, e in zip(q, el) if e])) for s, q in seqs.items()}
        valid = {s for s, c in cov.items() if c >= 0.80} | {TEMPLATE}
        frac = len(valid & set(placentals)) / len(placentals)
        if frac < 0.60:
            audit.append({'symbol': g, 'status': 'gene_filter', 'placental_fraction': frac})
            continue
        names = sorted(valid & set(species))
        tasks.append((g, names, placentals, asr.newick(asr.prune_tree(base_tree, {k: k for k in names}))))
        audit.append({'symbol': g, 'status': 'ok', 'placental_fraction': frac, 'n_taxa': len(names),
                      'eligible_sites': int(el.sum())})
    with ThreadPoolExecutor(max_workers=16) as pool:
        res = list(pool.map(run_asr, tasks))
    print('IQ-TREE', pd.Series([r['status'] for r in res]).value_counts().to_dict(), flush=True)

    model = load_model()
    rows, kept, post = [], [], []
    for r in res:
        if r['status'] != 'ok':
            continue
        g = r['symbol']
        probs = r['probs']
        anc = ''.join(AA20[i] for i in probs.argmax(axis=1))
        if not (ENERGY / f'{g}.npz').exists():
            d = parse(pdbs[g])
            if len(d['seq_chain_A']) != len(anc):
                continue
            np.savez_compressed(ENERGY / f'{g}.npz', e_anc=site_energies(model, d, anc), ancestor=anc)
        human_plddt = np.load(DATA / 'energies' / f'{g}.npz')['plddt']
        rows += gene_metrics(g, placentals, probs, pdbs[g], human_plddt)
        kept.append(g)
        post.append(r['mean_max_posterior'])
    met = pd.DataFrame(rows)
    met.to_csv(RESULTS / 'E6_species_gene_metrics.csv.gz', index=False)
    pd.DataFrame(audit).merge(pd.DataFrame([{k: v for k, v in r.items() if k != 'probs'} for r in res]),
                              on='symbol', how='left', suffixes=('', '_asr')).to_csv(RESULTS / 'E6_audit.csv', index=False)

    # Species x gene arrays for the new metric (two variants) and for E5 (human and platypus templates)
    sp = pd.read_csv(DATA / 'species_final.csv')
    sp = sp[sp.eligible & sp.placental].sort_values('ensembl_name').reset_index(drop=True)
    pl = sp.ensembl_name.tolist()
    prim = (sp.order == 'Primates').to_numpy()
    A = {}
    for v in ['main', 'human_sites']:
        for q in 'DK':
            A[q, v] = met.pivot(index='species', columns='symbol', values=f'{q}_{v}').reindex(index=pl, columns=kept).to_numpy()
    e5_arrays = {k: [] for k in ['human', 'platypus', 'K']}
    for g in kept:
        o, _ = e5.gene_sums(g, pl)
        for k in e5_arrays:
            e5_arrays[k].append(o[k])
    E5 = {k: np.column_stack(v) for k, v in e5_arrays.items()}

    def metrics(idx):
        m = {v: np.nansum(A['D', v][:, idx], 1) / np.nansum(A['K', v][:, idx], 1) for v in ['main', 'human_sites']}
        Kb = np.nansum(E5['K'][:, idx], 1)
        m.update({f'E5_{k}': np.nansum(E5[k][:, idx], 1) / Kb for k in ['human', 'platypus']})
        return m

    allidx = np.arange(len(kept))
    M = metrics(allidx)
    off = {k: float(offset(v, prim)) for k, v in M.items()}
    rng = np.random.default_rng(SEED_BOOT)
    boot = {k: [] for k in ['main', 'human_sites', 'E5_human', 'E5_platypus',
                            'main_minus_E5_platypus', 'main_minus_E5_human', 'human_sites_minus_main']}
    for _ in range(B_GENES):
        mb = metrics(rng.integers(0, len(kept), len(kept)))
        o = {k: offset(v, prim) for k, v in mb.items()}
        for k in ['main', 'human_sites', 'E5_human', 'E5_platypus']:
            boot[k].append(o[k])
        boot['main_minus_E5_platypus'].append(o['main'] - o['E5_platypus'])
        boot['main_minus_E5_human'].append(o['main'] - o['E5_human'])
        boot['human_sites_minus_main'].append(o['human_sites'] - o['main'])
    ci = {k: [float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))] for k, v in boot.items()}

    tree = Phylo.read(RAW / 'mcc_node.tre', 'nexus')
    model_p = PGLS(tree_covariance(tree, sp.tree_label.tolist()))
    y = np.log10(sp.max_lifespan_yr.to_numpy())
    lm = np.log10(sp.adult_mass_g.to_numpy())
    zf = lambda x: (x - x.mean()) / x.std(ddof=1)
    design = lambda m: np.column_stack([np.ones(len(y)), lm, zf(m)])
    fit, b = bootstrap_lrt(model_p, y, design(M['main']), 2, 999, SEED_P1)
    rank = lambda m: int(pd.Series(m).rank().iloc[pl.index('homo_sapiens')])
    zs = pd.Series(zf(M['main']), index=pl)
    out = {
        'amendment': 'E6 (post hoc, exploratory; written before computing)',
        'n_genes_input': len(genes), 'n_genes_retained': len(kept), 'n_species': len(pl),
        'audit': pd.DataFrame(audit).status.value_counts().to_dict(),
        'mean_max_posterior': float(np.mean(post)),
        'eligible_sites_total': int(sum(a['eligible_sites'] for a in audit if a['symbol'] in set(kept))),
        'primate_offset_SD': off, 'primate_offset_ci95': {k: ci[k] for k in ['main', 'human_sites', 'E5_human', 'E5_platypus']},
        'differences': {'main_minus_E5_platypus': {'estimate': off['main'] - off['E5_platypus'], 'ci95': ci['main_minus_E5_platypus']},
                        'main_minus_E5_human': {'estimate': off['main'] - off['E5_human'], 'ci95': ci['main_minus_E5_human']},
                        'human_sites_minus_main': {'estimate': off['human_sites'] - off['main'], 'ci95': ci['human_sites_minus_main']}},
        'human_rank_of_57': {k: rank(v) for k, v in M.items()},
        'human_z_main': float(zs['homo_sapiens']),
        'lowest_5_main': zs.sort_values().index[:5].tolist(),
        'by_order_mean_z_main': pd.Series(zf(M['main'])).groupby(sp.order.to_numpy()).mean().round(3).to_dict(),
        'r_with_E5': {k: float(np.corrcoef(M['main'], M[k])[0, 1]) for k in ['E5_human', 'E5_platypus', 'human_sites']},
        'P1_main': {**summarize(fit, 2), **b},
        'P1_wald': {k: summarize(model_p.fit(y, design(v)), 2) for k, v in M.items()},
    }
    pd.DataFrame({'ensembl_name': pl, 'order': sp.order, **{f'M_{k}': v for k, v in M.items()}}) \
        .to_csv(RESULTS / 'E6_species_metrics.csv', index=False)
    (RESULTS / 'E6_msa_platypus.json').write_text(json.dumps(out, indent=2, default=float))
    print(json.dumps(out, indent=1, default=float))


if __name__ == '__main__':
    main()
